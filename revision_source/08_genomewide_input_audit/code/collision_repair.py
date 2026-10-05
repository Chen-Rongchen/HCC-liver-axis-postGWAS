"""Order-independent collision handling and hash-bound frozen-cache overlays.

Only deterministic conflicting mapped keys are repaired here. Reverse-chain
ambiguity remains a separately labelled sensitivity in replay_scores.py.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from audit_inputs import BASE, ROOT, UPDATE, hashes

FIELDS=['p','beta','se','eaf','n']

def merge_mapped_records(arrays, destinations, values, effect_tokens):
    """Merge valid records, excluding all members of a conflicting target key.

    A conflict is sticky across later batches. Repeated indices are handled
    individually; fresh unique indices retain vectorized performance.
    """
    dst=np.asarray(destinations,dtype=np.int64);effect=np.asarray(effect_tokens)
    unique,count=np.unique(dst,return_counts=True)
    repeated=np.isin(dst,unique[count>1])|np.isfinite(arrays['p'][dst])|arrays['conflict'][dst].astype(bool)
    fresh=~repeated
    for field in FIELDS:arrays[field][dst[fresh]]=np.asarray(values[field])[fresh]
    arrays['effect_token'][dst[fresh]]=effect[fresh]
    for j in np.flatnonzero(repeated):
        k=dst[j]
        if arrays['conflict'][k]:continue
        if np.isfinite(arrays['p'][k]):
            same=all(arrays[f][k]==values[f][j] for f in FIELDS) and arrays['effect_token'][k]==effect[j]
            if not same:
                arrays['conflict'][k]=1
                for f in FIELDS:arrays[f][k]=np.nan
                arrays['effect_token'][k]=-1
        else:
            for f in FIELDS:arrays[f][k]=values[f][j]
            arrays['effect_token'][k]=effect[j]

def load_corrected_cache(accession, chromosome):
    """Return repaired in-memory arrays without modifying the original cache."""
    directory=BASE/'results'/accession
    provenance=pd.read_csv(directory/'source_hashes.tsv',sep='\t').set_index('path').sha256
    overlay=pd.read_csv(directory/'collision_repair_overlay.tsv',sep='\t')
    proof=pd.read_csv(directory/'cache_discrepancies.tsv',sep='\t')
    expected=proof.loc[proof.field.eq('p'),['chr','cache_index']].to_numpy(int)
    if not np.array_equal(overlay[['chr','cache_index']].to_numpy(int),expected):
        raise ValueError('Overlay no longer matches the full raw replay')
    out={}
    for f in FIELDS+['effect_token','conflict']:
        path=UPDATE/f'05_gene_level/adapted_variant_cache/{accession}/chr{chromosome:02d}_{f}.npy'
        if hashes(path)['sha256']!=provenance.loc[str(path.relative_to(ROOT))]:
            raise ValueError('Frozen cache changed; overlay is not applicable')
        out[f]=np.load(path).copy()
    idx=overlay.loc[overlay.chr.eq(chromosome),'cache_index'].to_numpy(int)
    for f in FIELDS:out[f][idx]=np.nan
    out['effect_token'][idx]=-1;out['conflict'][idx]=1
    return out
