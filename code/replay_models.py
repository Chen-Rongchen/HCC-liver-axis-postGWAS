"""Validate archived matrices by default; optionally reproduce EXISTING model fits."""
from pathlib import Path
import argparse, importlib.util, json
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('engine',ROOT/'code/model_engine.py')
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
p=argparse.ArgumentParser();p.add_argument('--refit',action='store_true',help='Reproduce saved coefficients/cluster-t CIs, not bootstrap or raw GWAS scoring')
args=p.parse_args();settings=json.loads((ROOT/'config/liver_only_model_settings.json').read_text())
expected=pd.read_csv(ROOT/'derived/model_estimates.tsv',sep='\t')
checks=[]
for branch in ['boundary','neutral']:
    for ancestry in ['EUR','EAS']:
        d=pd.read_csv(ROOT/f'derived/{branch}_{ancestry}_common_model_sample.tsv',sep='\t',float_precision='round_trip')
        assert len(d)==1652 and d.gene_symbol.is_unique
        for model in ['M0','M1','M2','M3']:
            x=engine.design(d,model,settings[branch])
            archived=pd.read_csv(ROOT/f'derived/{branch}_{ancestry}_{model}_design.tsv',sep='\t',float_precision='round_trip')
            assert np.allclose(x,archived[x.columns],rtol=1e-13,atol=1e-13)
            assert np.allclose(d.Y,archived.Y,rtol=0,atol=0)
            if args.refit:
                for scope in (['support','equal_block_weight'] if model=='M3' else ['support']):
                    groups=d[f'{ancestry}_ld_block'];weights=np.ones(len(d)) if scope=='support' else 1/d.groupby(f'{ancestry}_ld_block').gene_symbol.transform('size')
                    fit=sm.WLS(d.Y,x,weights=weights).fit(cov_type='cluster',cov_kwds={'groups':groups,'use_correction':True})
                    b=fit.params['high'];se=fit.bse['high'];q=t.ppf(.975,groups.nunique()-1)
                    e=expected.query('branch==@branch and ancestry==@ancestry and model==@model and scope==@scope').iloc[0]
                    assert np.allclose([b,se,b-q*se,b+q*se],[e.estimate,e.cluster_se,e.ci_low,e.ci_high],rtol=1e-7,atol=1e-8)
            checks.append(branch+'_'+ancestry+'_'+model)
print(json.dumps({'matrix_checks':len(checks),'refit_executed':args.refit,'raw_GWAS_pipeline_executed':False,'status':'PASS'},indent=2))
