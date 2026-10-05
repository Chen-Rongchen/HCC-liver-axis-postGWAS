"""Replay existing scores and quantify deterministic cache-collision repairs.

Original data remain read-only. The roundtrip exclusion is a separately labelled
conservative mapping sensitivity, not a claim to resolve physical identities.
"""
from pathlib import Path
from collections import defaultdict
import json, math
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm, fisher_exact, t
from audit_inputs import BASE, ROOT, MAP, UPDATE, hashes

RES=BASE/'results'; REL=ROOT/'github_release/HCC-liver-axis-postGWAS'
sources={}; checks=[]
def record(p):
    p=Path(p)
    if str(p) not in sources: sources[str(p)]=hashes(p)
    return p
def read(p): return pd.read_csv(record(p),sep='\t',float_precision='round_trip')
def arr(p): return np.load(record(p),mmap_mode='r')
def check(name, ok, **detail):
    checks.append(dict(name=name,passed=bool(ok),**detail))
    if not ok: raise AssertionError(checks[-1])

def score(p, boundary):
    p=np.asarray(p,float);p=p[np.isfinite(p)&(p>=0)&(p<=1)]
    if not len(p): return np.nan
    p=np.clip(p,1e-300,1);n=len(p)
    if n==1:return float(p[0])
    v=np.zeros(n);lo=p<.5;hi=(p>.5)&(p<1)
    v[lo]=1/np.tan(np.pi*p[lo]);v[hi]=-1/np.tan(np.pi*(1-p[hi]))
    if boundary:v[p==1]=0 if n==2 else -1/math.tan(math.pi/n)
    return max(math.atan2(1,math.fsum((v/n).tolist()))/math.pi,1e-300)

def score_tests():
    import mpmath as mp
    for boundary in [False,True]:
        for p in [[1.],[1.,1.],[1.,1.,1.],[1.,.001,.9],[1.,1e-300,.9],[.5,.5],[0.,1.,.2],[1.,np.nextafter(1.,0.),.01]]:
            with mp.workdps(400):
                v=[max(mp.mpf(float(x)),mp.mpf(1e-300)) for x in p];n=len(v)
                if n==1:ref=float(v[0])
                else:
                    terms=[mp.cot(mp.pi*(1-mp.mpf(1)/n)) if boundary and x==1 else mp.mpf(0) if x==1 else mp.cot(mp.pi*x) for x in v]
                    ref=max(float(mp.atan2(1,mp.fsum(terms)/n)/mp.pi),1e-300)
            check(f'ACAT_reference_{boundary}_{p}',math.isclose(score(p,boundary),ref,rel_tol=1e-10,abs_tol=0))

def ranks(df):
    for mode in ['boundary','neutral']:
        vals=df[mode+'_score'].to_numpy();valid=np.isfinite(vals);n=int(valid.sum())
        r=rankdata(vals[valid],method='average');fraction=np.full(len(df),np.nan);y=fraction.copy()
        fraction[valid]=r/n;y[valid]=norm.ppf((n-r+.5)/n)
        df[mode+'_rank_fraction']=fraction;df[mode+'_Y']=y
    return df

def main():
    score_tests();genes=read(MAP/'gene_registry.tsv');runs=np.zeros(len(genes),int)
    for ch in range(1,23):
        gi=arr(MAP/f'chr{ch:02d}_gene_index.npy');starts=np.r_[0,np.flatnonzero(np.diff(gi))+1]
        runs+=np.bincount(gi[starts],minlength=len(genes))
    summary=[];models=[];enrichment=[];membership=[];all_changes=[]
    masks=read(ROOT/'analysis_hbsn_restructure/35_seed_region_masking/out_03_membership_footprint/gene_mask_class_A_gene_body.tsv').set_index('gene').mask_class
    excluded=~masks.isin(['retained','retained_tss_unresolved'])
    for ancestry,accession in [('EUR','GCST90860790'),('EAS','GCST90860791')]:
        directory=RES/accession;old=read(REL/f'derived/{ancestry}_hcc_rank_scores.tsv').set_index('gene_symbol')
        delta=read(directory/'cache_discrepancies.tsv');bad=delta[delta.field.eq('p')]
        exceptions=read(directory/'mapping_and_numeric_exceptions.tsv.gz')
        conflicts=exceptions.merge(bad,left_on=['chromosome','cache_index'],right_on=['chr','cache_index'])
        check(ancestry+'_two_distinct_sources_per_conflict',conflicts.groupby(['chr','cache_index']).base_pair_location.nunique().eq(2).all())
        check(ancestry+'_conflict_P_values_differ',conflicts.groupby(['chr','cache_index']).p_value.nunique().eq(2).all())
        conflicts.to_csv(directory/'conflicting_source_records.tsv',sep='\t',index=False)
        # A compact overlay implements exclusion; no copied full input cache.
        overlay=bad[['chr','cache_index']].copy();overlay['action']='exclude_conflicting_mapped_key'
        overlay.to_csv(directory/'collision_repair_overlay.tsv',sep='\t',index=False)
        variants=['frozen','collision_repaired','roundtrip_exclusion_sensitivity']
        rows={v:{} for v in variants};parts={v:defaultdict(list) for v in variants}
        totals=np.zeros(len(genes),int)
        def add(v,g,p):
            finite=p[np.isfinite(p)]
            rows[v][g]=dict(gene_symbol=genes.iloc[g].gene,nsnps_matched=len(finite),nsnps_total=int(totals[g]),
                           boundary_score=score(finite,True),neutral_score=score(finite,False),n_P1=int((finite==1).sum()))
        for ch in range(1,23):
            gi=arr(MAP/f'chr{ch:02d}_gene_index.npy');vi=arr(MAP/f'chr{ch:02d}_variant_index.npy')
            original=arr(UPDATE/f'05_gene_level/adapted_variant_cache/{accession}/chr{ch:02d}_p.npy')
            corrected=np.asarray(original).copy();corrected[bad.loc[bad.chr.eq(ch),'cache_index'].to_numpy(int)]=np.nan
            strict=np.asarray(original).copy();strict[exceptions.loc[exceptions.chromosome.eq(ch),'cache_index'].to_numpy(int)]=np.nan
            starts=np.r_[0,np.flatnonzero(np.diff(gi))+1];stops=np.r_[starts[1:],len(gi)]
            for start,stop in zip(starts,stops):
                g=int(gi[start]);indices=vi[start:stop];totals[g]+=stop-start
                for v,p in zip(variants,[original,corrected,strict]):
                    values=p[indices]
                    if runs[g]==1:add(v,g,values)
                    else:parts[v][g].append(values)
            print(ancestry,'score chr',ch,flush=True)
        for v in variants:
            for g,p in parts[v].items():add(v,g,np.concatenate(p))
        check(ancestry+'_mapping_totals',np.array_equal(totals,genes.nsnps_cis))
        frames={v:ranks(pd.DataFrame(rows[v].values()).sort_values('gene_symbol').set_index('gene_symbol')) for v in variants}
        frozen=frames['frozen'].loc[old.index]
        for col in ['boundary_score','neutral_score','boundary_Y','neutral_Y','boundary_rank_fraction','neutral_rank_fraction','nsnps_matched','nsnps_total','n_P1']:
            check(ancestry+'_frozen_replay_'+col,np.allclose(frozen[col],old[col],rtol=1e-10,atol=0,equal_nan=True),max_abs_diff=float(np.nanmax(np.abs(frozen[col]-old[col]))))
        for v,frame in frames.items():
            frame.to_csv(directory/f'{v}_gene_scores.tsv.gz',sep='\t',float_format='%.17g')
            changed=~np.isclose(frame.boundary_score,frozen.loc[frame.index].boundary_score,rtol=1e-10,atol=0,equal_nan=True)
            summary.append(dict(ancestry=ancestry,version=v,evaluable=int(frame.boundary_score.notna().sum()),changed_boundary_scores=int(changed.sum()),
                changed_boundary_top5=int((frame.boundary_rank_fraction.le(.05)!=frozen.loc[frame.index].boundary_rank_fraction.le(.05)).sum())))
            for mode in ['boundary','neutral']:
                pop=read(REL/f'derived/{mode}_liver_only_population.tsv').set_index('gene_symbol')
                sample=read(REL/f'derived/{mode}_{ancestry}_common_model_sample.tsv').set_index('gene_symbol')
                admitted=pop.support_population & pop.index.to_series().map(frame[mode+'_Y']).notna() & pop[f'{ancestry}_ld_block'].notna()
                check(f'{ancestry}_{mode}_{v}_common_sample_members',set(pop.index[admitted])==set(sample.index))
                cand=pop.index[(~pop.seed_gene)&pop.B.ge(4)&pop.index.to_series().map(frame[mode+'_rank_fraction']).le(.05)]
                oldcand=pop.index[(~pop.seed_gene)&pop.B.ge(4)&pop.index.to_series().map(frozen[mode+'_rank_fraction']).le(.05)]
                membership.append(dict(ancestry=ancestry,mode=mode,version=v,model_members=len(sample),candidate_count=len(cand),candidate_added=len(set(cand)-set(oldcand)),candidate_removed=len(set(oldcand)-set(cand))))
                if ancestry=='EUR' and mode=='boundary' and v=='frozen':
                    events=read(REL/'derived/support_events.tsv')
                    check('original_283_candidate_set',set(cand)==set(events.gene_symbol.unique()) and len(cand)==283)
                if mode!='boundary':continue
                nonseed=pop[~pop.seed_gene & pop.index.to_series().map(frame.boundary_Y).notna()]
                for part,sub in [('original',nonseed),('retained',nonseed[~nonseed.index.to_series().map(excluded)])]:
                    nom=sub.liver_nomination_flag;top=sub.index.to_series().map(frame.boundary_rank_fraction).le(.05)
                    cells=[int((nom&top).sum()),int((nom&~top).sum()),int((~nom&top).sum()),int((~nom&~top).sum())]
                    odds,p=fisher_exact(np.array(cells).reshape(2,2))
                    enrichment.append(dict(ancestry=ancestry,version=v,part=part,nominee_top5=cells[0],nominee_total=sum(cells[:2]),background_top5=cells[2],background_total=sum(cells[2:]),OR=odds,fisher_p=p))
                design=read(REL/f'derived/boundary_{ancestry}_M2_design.tsv')
                design['Y']=design.gene_symbol.map(frame.boundary_Y)
                cols=['intercept','high','S','log_gene_length','log_cis_snp_count','log1p_adjacent_expression','stratum5']
                for part,sub in [('original',design),('retained',design[~design.gene_symbol.map(excluded)])]:
                    X=sub[cols].to_numpy(float);y=sub.Y.to_numpy();b=np.linalg.lstsq(X,y,rcond=None)[0]
                    groups=sub[f'{ancestry}_ld_block'].to_numpy();uniq=np.unique(groups);n,k=X.shape;G=len(uniq)
                    scores=X*(y-X@b)[:,None];meat=np.zeros((k,k))
                    for g in uniq:u=scores[groups==g].sum(axis=0);meat+=np.outer(u,u)
                    bread=np.linalg.inv(X.T@X);V=bread@meat@bread*G/(G-1)*(n-1)/(n-k);se=np.sqrt(V[1,1]);q=t.ppf(.975,G-1)
                    models.append(dict(ancestry=ancestry,version=v,part=part,n=n,blocks=G,estimate=b[1],cluster_se=se,ci_low=b[1]-q*se,ci_high=b[1]+q*se))
            delta_frame=frame.join(frozen[['boundary_score','boundary_Y','boundary_rank_fraction']],rsuffix='_frozen')
            delta_frame['version']=v;delta_frame['ancestry']=ancestry
            all_changes.append(delta_frame[changed])
    for name,values in [('score_impact',summary),('M2_impact',models),('enrichment_impact',enrichment),('membership_impact',membership)]:
        pd.DataFrame(values).to_csv(RES/f'{name}.tsv',sep='\t',index=False,float_format='%.17g')
    pd.concat(all_changes).to_csv(RES/'changed_gene_scores.tsv',sep='\t',float_format='%.17g')
    pd.DataFrame(sources.values()).to_csv(RES/'score_replay_sources.tsv',sep='\t',index=False)
    (RES/'score_replay_checks.json').write_text(json.dumps(dict(checks=checks,all_passed=all(x['passed'] for x in checks),new_model_specifications=0,existing_M2_refits=len(models)),indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
