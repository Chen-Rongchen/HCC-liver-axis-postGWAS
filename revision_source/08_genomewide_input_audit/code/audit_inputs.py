"""Independent whole-genome raw-to-cache replay; original files are read-only."""
from pathlib import Path
from collections import defaultdict, Counter
import argparse,gzip,hashlib,json,math
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1];ROOT=BASE.parents[1];RES=BASE/'results'
UPDATE=ROOT/'analysis_update_2026_hcc_gwas';MAP=UPDATE/'03_qc/frozen_map_cache'
FIELDS={'p':'p_value','beta':'beta','se':'standard_error','eaf':'effect_allele_frequency','n':'n'}
EXPECTED_MD5={'GCST90860790':'13b3a1f9a4120747bbe55a7905439bd1','GCST90860791':'e0d4d1266d92271be54cd8bcebb78055'}

def hashes(p):
    a=hashlib.md5();b=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):a.update(chunk);b.update(chunk)
    return dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,md5=a.hexdigest(),sha256=b.hexdigest())

def segments(blocks,back=False):
    """Sweep intervals; unlike last-start lookup, detect nested/overlapping chains."""
    events=defaultdict(list)
    for idx,b in enumerate(blocks):
        start,stop=(b[2],b[2]+b[1]-b[0]) if back else b[:2]
        events[start].append((idx,1));events[stop].append((idx,-1))
    active=set();starts=[];counts=[];offsets=[]
    for pos in sorted(events):
        for idx,sign in events[pos]:
            if sign==1:active.add(idx)
            else:active.remove(idx)
        starts.append(pos);counts.append(len(active))
        if len(active)==1:
            b=blocks[next(iter(active))];offsets.append(b[0]-b[2] if back else b[2]-b[0])
        else:offsets.append(0)
    return np.array(starts,np.int64),np.array(counts,np.int32),np.array(offsets,np.int64)

def read_chains():
    out=defaultdict(list)
    with gzip.open(ROOT/'data/reference/liftover/hg38ToHg19.over.chain.gz','rt') as f:
        active=False
        for line in f:
            x=line.split()
            if not x:continue
            if x[0]=='chain':
                active=x[2]==x[7] and x[4]==x[9]=='+' and x[2][3:] in {str(i) for i in range(1,23)}
                if active:ch=int(x[2][3:]);t=int(x[5]);q=int(x[10])
            elif active:
                size=int(x[0]);out[ch].append((t,t+size,q))
                if len(x)==3:t+=size+int(x[1]);q+=size+int(x[2])
                else:active=False
    return {c:dict(blocks=sorted(b),forward=segments(b),reverse=segments(b,True)) for c,b in out.items()}

def query(seg,positions):
    starts,counts,offsets=seg
    loc=np.searchsorted(starts,positions,side='right')-1;safe=np.maximum(loc,0)
    n=np.where(loc>=0,counts[safe],0)
    return n,positions+offsets[safe]

def audit(accession):
    lock=json.loads((BASE/'protocol_lock.json').read_text())
    assert hashes(ROOT/lock['protocol'])['sha256']==lock['sha256']
    raw=UPDATE/f'02_inputs/raw_gwas/{accession}.tsv.gz';rawhash=hashes(raw)
    assert rawhash['md5']==EXPECTED_MD5[accession]
    directory=RES/accession;directory.mkdir(exist_ok=True)
    if (directory/'summary.json').exists() or (directory/'mapping_and_numeric_exceptions.tsv.gz').exists():
        raise FileExistsError('Audit output already exists; use a fresh audit directory to preserve the source trace and prevent duplicate exception logs')
    original=UPDATE/f'05_gene_level/adapted_variant_cache/{accession}'
    chains=read_chains();data={};provenance=[rawhash,hashes(ROOT/'data/reference/liftover/hg38ToHg19.over.chain.gz')]
    for ch in range(1,23):
        paths={f:original/f'chr{ch:02d}_{f}.npy' for f in [*FIELDS,'effect_token','conflict']}
        keypath=MAP/f'chr{ch:02d}_variant_keys.npy';tokenpath=MAP/f'chr{ch:02d}_allele_tokens.tsv'
        keys=np.load(keypath,mmap_mode='r');tokens=pd.read_csv(tokenpath,sep='\t',keep_default_na=False).set_index('allele').allele_code.to_dict()
        data[ch]=dict(keys=keys,tokens=tokens,old={f:np.load(p,mmap_mode='r') for f,p in paths.items()},
            replay={f:np.full(len(keys),np.nan) for f in FIELDS},effect=np.full(len(keys),-1,np.int32),
            conflict=np.zeros(len(keys),bool),invalid_seen=np.zeros(len(keys),bool),seen=np.zeros(len(keys),bool))
        provenance.extend(hashes(p) for p in [keypath,tokenpath,*paths.values()])
    counters=Counter();anomalies=[];lastcoord=-1
    # Save only exceptions, never another complete GWAS or cache copy.
    exception=directory/'mapping_and_numeric_exceptions.tsv.gz';header=True
    def process(d):
        nonlocal lastcoord,header
        if d.empty:return
        chrs=d.chromosome.to_numpy(int);pos=d.base_pair_location.to_numpy(np.int64)
        coord=chrs.astype(np.int64)*1000000000+pos
        assert coord[0]>=lastcoord and (np.diff(coord)>=0).all(),'Unsorted input: duplicate audit cannot be assumed complete'
        lastcoord=coord[-1];counters['raw_rows']+=len(d)
        ea=d.effect_allele.str.upper();oa=d.other_allele.str.upper()
        d=d.copy();d['lo']=np.minimum(ea.to_numpy(),oa.to_numpy());d['hi']=np.maximum(ea.to_numpy(),oa.to_numpy())
        duplicate=d.duplicated(['chromosome','base_pair_location','lo','hi'],keep=False)
        if duplicate.any():
            keep=pd.Series(True,index=d.index)
            for _,g in d[duplicate].groupby(['chromosome','base_pair_location','lo','hi']):
                identical=len(g.drop_duplicates())==1
                counters['source_duplicate_groups']+=1
                if identical:keep.loc[g.index[1:]]=False;counters['source_identical_duplicate_rows_removed']+=len(g)-1
                else:keep.loc[g.index]=False;counters['source_conflicting_duplicate_rows_removed']+=len(g)
            d=d[keep];ea=d.effect_allele.str.upper();oa=d.other_allele.str.upper()
        p=d.p_value.to_numpy(dtype=float)
        valid=(np.isfinite(p)&(p>=0)&(p<=1)&np.isfinite(d.beta)&np.isfinite(d.standard_error)&d.standard_error.gt(0)
               &np.isfinite(d.effect_allele_frequency)&d.effect_allele_frequency.between(0,1)&np.isfinite(d.n)&d.n.gt(0))
        d['p_numeric']=p;d['numeric_valid']=valid;d['ea']=ea;d['oa']=oa
        counters['numeric_invalid_raw_rows']+=int((~valid).sum());counters['P_zero_or_underflow_rows']+=int((p==0).sum())
        for ch,g in d.groupby('chromosome',sort=False):
            ch=int(ch)
            if ch not in data:counters['nonautosomal_rows']+=len(g);continue
            x=g.base_pair_location.to_numpy(np.int64)-1;b=chains[ch]
            nf,safe37=query(b['forward'],x)
            # Replay historical mapping algebra independently, while auditing full uniqueness.
            blocks=np.asarray(b['blocks'],np.int64);ix=np.searchsorted(blocks[:,0],x,side='right')-1;loc=np.maximum(ix,0)
            mapped=(ix>=0)&(x<blocks[loc,1]);p37=x-blocks[loc,0]+blocks[loc,2]+1
            counters['same_chrom_plus_mapped_rows']+=int(mapped.sum())
            counters['forward_multiple_raw_rows']+=int((nf>1).sum())
            counters['legacy_missed_unique_mapping_rows']+=int(((nf==1)&~mapped).sum())
            counters['legacy_unique_mapping_disagreements']+=int(((nf==1)&mapped&(p37!=safe37+1)).sum())
            nb,back38=query(b['reverse'],p37-1)
            dat=data[ch];tokens=dat['tokens'];keys=dat['keys']
            a=g.ea.map(tokens).fillna(-1).to_numpy(np.int64);z=g.oa.map(tokens).fillna(-1).to_numpy(np.int64)
            candidate=mapped&(p37>0)&(a>=0)&(z>=0)&(a!=z)
            encoded=(np.maximum(p37,0).astype(np.uint64)<<np.uint64(34))|(np.minimum(np.maximum(a,0),np.maximum(z,0)).astype(np.uint64)<<np.uint64(17))|np.maximum(np.maximum(a,0),np.maximum(z,0)).astype(np.uint64)
            loc=np.searchsorted(keys,encoded);sloc=np.minimum(loc,len(keys)-1)
            match=candidate&(loc<len(keys))&(keys[sloc]==encoded)
            good=match&g.numeric_valid.to_numpy();bad=match&~g.numeric_valid.to_numpy()
            counters['exact_frozen_key_raw_matches']+=int(match.sum());counters['valid_frozen_key_raw_matches']+=int(good.sum())
            counters['invalid_frozen_key_raw_matches']+=int(bad.sum())
            ambiguous=good&((nf!=1)|(nb!=1)|(back38!=x))
            counters['accepted_ambiguous_mapping_rows']+=int(ambiguous.sum())
            dat['invalid_seen'][sloc[bad]]=True
            interesting=bad|ambiguous
            if interesting.any():
                out=g.loc[interesting,['chromosome','base_pair_location','rsid','effect_allele','other_allele','p_value','beta','standard_error','effect_allele_frequency','n']].copy()
                out['pos37']=p37[interesting];out['forward_hits']=nf[interesting];out['reverse_hits']=nb[interesting]
                out['numeric_invalid']=bad[interesting];out['accepted_ambiguous']=ambiguous[interesting];out['cache_index']=sloc[interesting]
                with gzip.open(exception,'at') as f:out.to_csv(f,sep='\t',index=False,header=header)
                header=False
            src=np.flatnonzero(good);dst=sloc[good]
            vals={f:g[c].to_numpy(float)[good] for f,c in FIELDS.items()};vals['p']=np.where(vals['p']==0,1e-300,vals['p'])
            # Most rows are unique. Process repeated target indices individually to
            # avoid numpy assignment hiding within-chunk or inter-chunk collisions.
            unique,ct=np.unique(dst,return_counts=True);repeated=np.isin(dst,unique[ct>1])|dat['seen'][dst]
            fresh=~repeated
            for f in FIELDS:dat['replay'][f][dst[fresh]]=vals[f][fresh]
            dat['effect'][dst[fresh]]=a[src[fresh]];dat['seen'][dst[fresh]]=True
            for j in np.flatnonzero(repeated):
                k=dst[j]
                if dat['seen'][k]:
                    same=all(dat['replay'][f][k]==vals[f][j] for f in FIELDS) and dat['effect'][k]==a[src[j]]
                    dat['conflict'][k]|=not same;counters['repeated_mapped_keys']+=1
                else:
                    for f in FIELDS:dat['replay'][f][k]=vals[f][j]
                    dat['effect'][k]=a[src[j]];dat['seen'][k]=True
    pending=None
    for chunk in pd.read_csv(raw,sep='\t',chunksize=300000,float_precision='round_trip',dtype={'p_value':str,'effect_allele':str,'other_allele':str}):
        d=pd.concat([pending,chunk],ignore_index=True) if pending is not None else chunk
        end=d.iloc[-1];mask=d.chromosome.eq(end.chromosome)&d.base_pair_location.eq(end.base_pair_location)
        pending=d[mask].copy();process(d[~mask])
        print(accession,'scanned',counters['raw_rows'],flush=True)
    process(pending)
    rows=[];differences=[];impact_indices={}
    for ch,dat in data.items():
        conflict=dat['conflict']
        for f in FIELDS:dat['replay'][f][conflict]=np.nan
        dat['effect'][conflict]=-1
        old=dat['old'];expected=dat['replay'];impact=np.zeros(len(dat['keys']),bool)
        for f in FIELDS:
            av=old[f];bv=expected[f];missing=np.isfinite(av)!=np.isfinite(bv)
            finite=np.isfinite(av)&np.isfinite(bv)
            close=np.isclose(av,bv,rtol=1e-12,atol=0 if f=='p' else 1e-15,equal_nan=True)
            numeric=finite&~close;diff=missing|numeric
            if f=='p':impact|=diff
            rel=np.abs(av[finite]-bv[finite])/np.maximum(np.abs(bv[finite]),1e-300)
            rows.append(dict(chr=ch,field=f,keys=len(av),expected_finite=int(np.isfinite(bv).sum()),cache_finite=int(np.isfinite(av).sum()),missing_mismatches=int(missing.sum()),numeric_mismatches=int(numeric.sum()),max_relative_error=float(rel.max()) if len(rel) else 0))
            for k in np.flatnonzero(diff):differences.append(dict(chr=ch,cache_index=k,field=f,cache_value=float(av[k]),raw_replay_value=float(bv[k])))
        for name,new in [('effect_token',dat['effect']),('conflict',conflict)]:
            diff=old[name]!=new
            rows.append(dict(chr=ch,field=name,keys=len(new),expected_finite=len(new),cache_finite=len(new),missing_mismatches=0,numeric_mismatches=int(diff.sum()),max_relative_error=0))
        impact_indices[ch]=np.flatnonzero(impact)
        counters['invalid_only_unrepresented_keys']+=int((dat['invalid_seen']&~dat['seen']).sum())
        counters['unique_valid_keys']+=int(np.isfinite(expected['p']).sum())
    comparison=pd.DataFrame(rows);comparison.to_csv(directory/'cache_comparison.tsv',sep='\t',index=False)
    pd.DataFrame(differences,columns=['chr','cache_index','field','cache_value','raw_replay_value']).to_csv(directory/'cache_discrepancies.tsv',sep='\t',index=False)
    genes=pd.read_csv(MAP/'gene_registry.tsv',sep='\t');affected=np.zeros(len(genes),bool)
    for ch,idx in impact_indices.items():
        if not len(idx):continue
        gi=np.load(MAP/f'chr{ch:02d}_gene_index.npy',mmap_mode='r');vi=np.load(MAP/f'chr{ch:02d}_variant_index.npy',mmap_mode='r')
        affected[np.unique(gi[np.isin(vi,idx)])]=True
    genes.loc[affected].to_csv(directory/'P_cache_affected_genes.tsv',sep='\t',index=False)
    total=int(comparison[['missing_mismatches','numeric_mismatches']].to_numpy().sum())
    result=dict(accession=accession,status='pass_frozen_P_cache_replay' if total==0 else 'deterministic_discrepancy_found',
        counters=dict(counters),cache_field_discrepancies=total,P_affected_genes=int(affected.sum()),
        mapping_uniqueness_pass=counters['accepted_ambiguous_mapping_rows']==0,
        raw_identity=rawhash,scope='raw file to all frozen cache keys; original numerical filter and P floor retained',new_models_run=0)
    (directory/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(provenance).to_csv(directory/'source_hashes.tsv',sep='\t',index=False)
    print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('accession',choices=list(EXPECTED_MD5));args=p.parse_args();audit(args.accession)
