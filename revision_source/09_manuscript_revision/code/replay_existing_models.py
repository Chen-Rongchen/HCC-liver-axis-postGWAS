"""Synchronize existing reported model outputs after deterministic input repair."""
from pathlib import Path
import importlib.util,json,shutil,hashlib
import pandas as pd
import numpy as np
BASE=Path(__file__).resolve().parents[1];ROOT=BASE.parents[1]
REL=ROOT/'github_release/HCC-liver-axis-postGWAS';OUT=BASE/'results/models';OUT.mkdir(parents=True,exist_ok=True)
source=REL/'code/model_engine.py'
spec=importlib.util.spec_from_file_location('existing_model_engine',source);engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
engine.OUT=OUT;engine.ROOT=ROOT;engine.BASE=ROOT/'analysis_hbsn_restructure'
for name in ['boundary_liver_only_population.tsv','neutral_liver_only_population.tsv']:
    shutil.copy2(REL/'derived'/name,OUT/name)
shutil.copy2(REL/'config/liver_only_model_settings.json',OUT/'liver_only_model_settings.json')
for anc,acc in [('EUR','GCST90860790'),('EAS','GCST90860791')]:
    old=pd.read_csv(REL/f'derived/{anc}_hcc_rank_scores.tsv',sep='\t',float_precision='round_trip').set_index('gene_symbol')
    fixed=pd.read_csv(ROOT/f'analysis_signal_resolution_c/08_genomewide_input_audit/results/{acc}/collision_repaired_gene_scores.tsv.gz',sep='\t',float_precision='round_trip').set_index('gene_symbol')
    for c in fixed.columns:old[c]=fixed.loc[old.index,c]
    old['mapping_fraction']=old.nsnps_matched/old.nsnps_total
    old=old.rename(columns={'min_p':'historical_min_p','max_n_author_weight':'historical_max_n_author_weight'})
    old['raw_n_policy']='Source METAL Weight is provenance only; beta scale and RSS N semantics are not certified'
    old.to_csv(OUT/f'{anc}_hcc_rank_scores.tsv',sep='\t',float_format='%.17g')
engine.analyze()
# Exact matrices and source populations are supplied separately from the journal files.
settings=json.loads((OUT/'liver_only_model_settings.json').read_text())
for branch in ['boundary','neutral']:
    for anc in ['EUR','EAS']:
        sample=engine.read(OUT/f'{branch}_{anc}_common_model_sample.tsv')
        for m in ['M0','M1','M2','M3']:
            x=engine.design(sample,m,settings[branch]);x.insert(0,'gene_symbol',sample.gene_symbol);x['Y']=sample.Y;x[f'{anc}_ld_block']=sample[f'{anc}_ld_block']
            engine.write(x,f'{branch}_{anc}_{m}_design.tsv')
previous=ROOT/'submission/EJHG'
if not previous.exists():previous=BASE/'baseline_package_files'
new=engine.read(OUT/'model_estimates.tsv');old=engine.read(previous/'Supplementary_Source_TSV/S09_model_estimates.tsv')
keys=['branch','ancestry','scope','model'];merged=old.merge(new,on=keys,suffixes=('_old','_new'),validate='one_to_one')
assert len(merged)==56 and (merged.n_old==merged.n_new).all()
merged.to_csv(OUT/'model_change_trace.tsv',sep='\t',index=False,float_format='%.17g')
lobo=[]
for branch in ['boundary','neutral']:
    for anc in ['EUR','EAS']:
        d=engine.read(OUT/f'lobo_{branch}_{anc}_support_M3.tsv');r=d.loc[d.delta.abs().idxmax()]
        lobo.append(dict(branch=branch,ancestry=anc,minimum_leave_out_estimate=d.estimate.min(),maximum_leave_out_estimate=d.estimate.max(),largest_influence_block=r.block,estimate_without_largest_influence=r.estimate,all_valid=bool(d.valid.all())))
engine.write(pd.DataFrame(lobo),'block_influence_summary.tsv')
(OUT/'replay_provenance.json').write_text(json.dumps(dict(original_engine=str(source.relative_to(ROOT)),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),models=56,new_model_specifications=0,all_engine_checks_passed=all(c['passed'] for c in engine.CHECKS),checks=len(engine.CHECKS)),indent=2)+'\n')
