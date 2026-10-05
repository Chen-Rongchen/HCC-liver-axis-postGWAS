"""Actual legacy adapter regression on all 49 real conflicting target keys."""
from collections import Counter
import importlib.util, json
import numpy as np
import pandas as pd
from audit_inputs import BASE, ROOT, MAP, hashes
from collision_repair import FIELDS,merge_mapped_records,load_corrected_cache

path=ROOT/'analysis_update_2026_hcc_gwas/code/qc_primary_2026_endpoints.py'
spec=importlib.util.spec_from_file_location('frozen_adapter',path);legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
chain=legacy.SameChromosomeChainMap(legacy.CHAIN)
def empty(n):
    return {**{f:np.full(n,np.nan) for f in FIELDS},'effect_token':np.full(n,-1,int),'conflict':np.zeros(n,np.uint8)}
checks=[]
def check(name,ok):
    checks.append(dict(name=name,passed=bool(ok)))
    if not ok:raise AssertionError(name)

for accession in ['GCST90860790','GCST90860791']:
    directory=BASE/'results'/accession
    data=pd.read_csv(directory/'conflicting_source_records.tsv',sep='\t',float_precision='round_trip')
    for (ch,k),g in data.groupby(['chr','cache_index']):
        # Instantiate without __init__: all writes are to one-element RAM arrays.
        scan=object.__new__(legacy.EndpointScan);scan.accession=accession;scan.stats=Counter();scan.duplicate_rows=[]
        scan.chain=chain
        keys=np.load(MAP/f'chr{ch:02d}_variant_keys.npy',mmap_mode='r')
        tokens=pd.read_csv(MAP/f'chr{ch:02d}_allele_tokens.tsv',sep='\t',keep_default_na=False).set_index('allele').allele_code.to_dict()
        scan.keys={ch:keys[[k]]};scan.positions={ch:np.array([g.pos37.iloc[0]])};scan.token_lookup={ch:tokens};scan.arrays={ch:empty(1)}
        scan.adapt(g.copy())
        check(f'{accession}_{ch}_{k}_legacy_overwrite',scan.arrays[ch]['conflict'][0]==0 and scan.arrays[ch]['p'][0]==g.cache_value.iloc[0])
        vals={f:g[c].to_numpy() for f,c in zip(FIELDS,['p_value','beta','standard_error','effect_allele_frequency','n'])};effect=g.effect_allele.map(tokens).to_numpy()
        for order in [np.array([0,1]),np.array([1,0])]:
            for split in [False,True]:
                a=empty(1)
                for ix in ([order[:1],order[1:]] if split else [order]):
                    merge_mapped_records(a,np.zeros(len(ix),int),{f:v[ix] for f,v in vals.items()},effect[ix])
                check(f'{accession}_{ch}_{k}_fixed_{order.tolist()}_{split}',a['conflict'][0]==1 and all(np.isnan(a[f][0]) for f in FIELDS) and a['effect_token'][0]==-1)
    # Exercise actual overlay reader on every chromosome with a correction.
    for ch in sorted(data.chr.unique()):
        patched=load_corrected_cache(accession,int(ch));indices=data.loc[data.chr.eq(ch),'cache_index'].unique()
        check(f'{accession}_{ch}_overlay',patched['conflict'][indices].all() and np.isnan(patched['p'][indices]).all())

# Identical duplicates stay valid; a conflict cannot be resurrected later.
a=empty(2);v={f:np.array([.25,.25,.7]) for f in FIELDS};v['se'][:]=1;v['n'][:]=100
merge_mapped_records(a,[0,0,1],v,[1,1,1])
check('identical_duplicates_retained',a['p'][0]==.25 and not a['conflict'].any())
merge_mapped_records(a,[0],{f:np.array([2.]) for f in FIELDS},[1])
merge_mapped_records(a,[0],{f:np.array([.25]) for f in FIELDS},[1])
check('conflict_sticky',a['conflict'][0]==1 and np.isnan(a['p'][0]))
(BASE/'results/collision_tests.json').write_text(json.dumps(dict(checks=checks,all_passed=True,legacy_source=hashes(path)),indent=2)+'\n')
print('PASS',len(checks),'collision checks')
