"""Finite post-audit reanalysis. All writes are confined to this directory."""
from pathlib import Path
import sys, json, hashlib, importlib.util, argparse, math
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm, t, fisher_exact
import statsmodels.api as sm

OUT=Path(__file__).resolve().parent; BASE=OUT.parent; ROOT=BASE.parent
TRAITS=['CIRR','LIVER_FAT','NAFLD','ALT','AST','GGT','ALP','BILI','PBC','PSC','AIH']
XCOL=['S','log_gene_length','log_cis_snp_count','log1p_adjacent_expression']
CHECKS=[]
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path); mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod;spec.loader.exec_module(mod);return mod
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for s in iter(lambda:f.read(8*1024*1024),b''):h.update(s)
    return h.hexdigest()
def read(path):return pd.read_csv(path,sep='\t',float_precision='round_trip')
def write(df,name):df.to_csv(OUT/name,sep='\t',index=False,na_rep='',float_format='%.17g')
def check(name,condition,detail=''):
    CHECKS.append(dict(check=name,passed=bool(condition),detail=str(detail)))
    if not condition:raise AssertionError(f'{name}: {detail}')
def rcs(z,knots):
    pos=lambda a:np.maximum(a,0)**3
    a,b=knots[-2:]
    return np.column_stack([z]+[(pos(z-k)-pos(z-a)*(b-k)/(b-a)+pos(z-b)*(a-k)/(b-a))/(knots[-1]-knots[0])**2 for k in knots[:-2]])

def prepare():
    """Only liver-side values and block metadata; do not load HCC statistics."""
    fp=load(BASE/'code/01_audit_identifiability.py','final_feature_paths')
    manifest=[];pop=[];ranges=[];settings={}
    for branch,path in [('boundary',BASE/'25_p1_boundary_sensitivity/boundary_liver_gene_registry.tsv'),('neutral',BASE/'23_liver_acat_audit/stable_liver_gene_registry.tsv')]:
        g=read(path);manifest.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
        check(branch+'_unique_universe',len(g)==18321 and g.gene_symbol.is_unique and g.seed_gene.sum()==56)
        f=g[['gene_symbol','chr','seed_gene','liver_nomination_flag','n_liver_traits_supported','n_domains_supported']].copy()
        f['B']=f.n_liver_traits_supported;f['high']=f.B.ge(4).astype(int)
        f['S']=1-g[[c+'_rank_pct' for c in TRAITS]].min(axis=1)
        for path,cols in [(fp.GENE_LENGTH,['gene','gene_length_bp']),(fp.CIS_SNP_COUNT,['gene','nsnps_cis']),
                         (fp.EXPRESSION,['gene_symbol','HCCDB18_adjacent_median']),
                         (fp.LD_BLOCK,['gene','EUR_ld_block','EAS_ld_block'])]:
            v=pd.read_csv(path,sep='\t',usecols=cols).rename(columns={'gene':'gene_symbol'})
            f=f.merge(v,on='gene_symbol',how='left',validate='one_to_one')
            manifest.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
        f['S_quintile']=np.nan
        nominees=f.liver_nomination_flag&~f.seed_gene
        f.loc[nominees,'S_quintile']=pd.qcut(f.loc[nominees,'S'],5,labels=[1,2,3,4,5]).astype(int).to_numpy()
        f['complete_X']=f.gene_length_bp.gt(0)&f.nsnps_cis.ge(0)&f.HCCDB18_adjacent_median.ge(0)
        f['log_gene_length']=np.log(f.gene_length_bp.where(f.gene_length_bp.gt(0)))
        f['log_cis_snp_count']=np.log1p(f.nsnps_cis)
        f['log1p_adjacent_expression']=np.log1p(f.HCCDB18_adjacent_median.where(f.HCCDB18_adjacent_median.ge(0)))
        eligible=nominees&f.complete_X&f.S_quintile.isin([4,5]);f['support_population']=False
        for q in [4,5]:
            d=f.loc[eligible&f.S_quintile.eq(q)].copy(); keep=pd.Series(True,index=d.index)
            for c in XCOL:
                a=d.loc[d.high.eq(0),c];b=d.loc[d.high.eq(1),c]
                lo=max(a.min(),b.min());hi=min(a.max(),b.max())
                ranges.append(dict(branch=branch,S_quintile=q,variable=c,lower=lo,upper=hi))
                keep &= d[c].between(lo,hi)
            f.loc[d.index[keep],'support_population']=True
            for phase,mask in [('before',pd.Series(True,index=d.index)),('retained',keep)]:
                for high in [0,1]:
                    z=d.loc[mask&d.high.eq(high)];pop.append(dict(branch=branch,S_quintile=q,phase=phase,high=high,n=len(z),EUR_blocks=z.EUR_ld_block.nunique()))
        d=f[f.support_population]
        settings[branch]={}
        for c in XCOL:
            settings[branch][c]={'mean':float(d[c].mean()),'sd':float(d[c].std(ddof=0))}
        for c in ['S','log_cis_snp_count']:
            z=(d[c]-settings[branch][c]['mean'])/settings[branch][c]['sd']
            knots=np.quantile(z,[.05,.35,.65,.95]);check(branch+c+'_knots',np.all(np.diff(knots)>0))
            settings[branch][c]['knots']=knots.tolist()
        write(f,f'{branch}_liver_only_population.tsv')
    write(pd.DataFrame(manifest).drop_duplicates(),'input_manifest.tsv');write(pd.DataFrame(pop),'population_before_outcome.tsv');write(pd.DataFrame(ranges),'marginal_support_ranges.tsv')
    (OUT/'liver_only_model_settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    frozen=['charter.md','liver_only_model_settings.json','boundary_liver_only_population.tsv','neutral_liver_only_population.tsv']
    (OUT/'freeze_manifest.json').write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
          HCC_loaded=False,scope='post_result_amendment_not_prospective_registration',hashes={n:sha(OUT/n) for n in frozen}),indent=2)+'\n')
    write(pd.DataFrame(CHECKS),'prepare_validation.tsv');print(pd.DataFrame(pop).to_string(index=False),flush=True)

def scores():
    frozen=json.loads((OUT/'freeze_manifest.json').read_text())
    for name,h in frozen['hashes'].items():check(name+'_frozen',sha(OUT/name)==h)
    mod=load(BASE/'25_p1_boundary_sensitivity/run_traits.py','final_boundary_functions')
    engine=load(ROOT/'analysis_update_2026_hcc_gwas/code/compute_primary_gene_level_acat.py','final_hcc_engine')
    genes=read(engine.MAP/'gene_registry.tsv');stats=[];refs=[]
    def score(values):
        a=mod.boundary(values);b=mod.neutral(values)
        if (values==1).any() and len(refs)<6:
            reference=mod.reference(values);check('mp400_actual_HCC',np.isclose(a,reference,rtol=1e-10,atol=1e-15))
            refs.append(dict(d=len(values),score=a,mp400=reference))
        return dict(boundary_score=a,neutral_score=b,n_P1=int((values==1).sum()))
    engine.acat=score
    for ancestry,accession in [('EUR','GCST90860790'),('EAS','GCST90860791')]:
        refs.clear();g=engine.compute_endpoint(accession,engine.ENDPOINTS[accession],genes)
        g=pd.concat([g.drop(columns='acat_p'),pd.DataFrame(g.acat_p.tolist())],axis=1).rename(columns={'gene':'gene_symbol'})
        old=read(BASE/f'22_acat_numerical_audit/{ancestry}_rank_score_deltas.tsv').set_index('gene_symbol')
        check(ancestry+'_neutral_replay',np.allclose(g.neutral_score,old.loc[g.gene_symbol,'stable_p'],equal_nan=True,rtol=1e-12,atol=0))
        check(ancestry+'_P1_absent_identical',np.allclose(g.loc[g.n_P1.eq(0),'boundary_score'],g.loc[g.n_P1.eq(0),'neutral_score'],equal_nan=True,rtol=1e-12,atol=0))
        for branch in ['boundary','neutral']:
            valid=g[f'{branch}_score'].notna();n=valid.sum();r=rankdata(g.loc[valid,f'{branch}_score'],method='average')
            g.loc[valid,f'{branch}_rank_fraction']=r/n
            g.loc[valid,f'{branch}_Y']=norm.ppf((n-r+.5)/n)
        stats.append(dict(ancestry=ancestry,evaluable=int(g.evaluable.sum()),P1_genes=int(g.n_P1.gt(0).sum()),
                          rank_changed=int(g.boundary_rank_fraction.ne(g.neutral_rank_fraction).where(g.evaluable,False).sum()),
                          top5_entered=int((g.boundary_rank_fraction.le(.05)&~g.neutral_rank_fraction.le(.05)).sum()),
                          top5_exited=int((~g.boundary_rank_fraction.le(.05)&g.neutral_rank_fraction.le(.05)).sum())))
        write(g,f'{ancestry}_hcc_rank_scores.tsv');write(pd.DataFrame(refs),f'{ancestry}_mp400_validation.tsv')
    write(pd.DataFrame(stats),'hcc_scoring_summary.tsv');write(pd.DataFrame(CHECKS),'score_validation.tsv')

def design(d,model,settings,exposure='high'):
    x=pd.DataFrame({'intercept':np.ones(len(d)),exposure:d[exposure].astype(float).to_numpy()},index=d.index)
    if model=='M0':return x
    cols=['S'] if model=='M1' else XCOL
    for c in cols:
        z=(d[c].to_numpy()-settings[c]['mean'])/settings[c]['sd']
        if model=='M3' and c in ['S','log_cis_snp_count']:
            basis=rcs(z,np.array(settings[c]['knots']))
            for j in range(basis.shape[1]):x[f'{c}_{j}']=basis[:,j]
        else:x[c]=z
    if d.S_quintile.nunique()>1:x['stratum5']=d.S_quintile.eq(5).astype(float).to_numpy()
    return x

def fit(d,x,ancestry,branch,scope,model,weights=None,bootstrap=False,lobo=False):
    y=d.Y.to_numpy();a=x.to_numpy(float);block=d[f'{ancestry}_ld_block'].astype(str).to_numpy()
    w=np.ones(len(d)) if weights is None else weights;G=np.unique(block);ng=len(G)
    check(f'{branch}_{ancestry}_{scope}_{model}_rank',np.linalg.matrix_rank(a)==a.shape[1])
    result=sm.WLS(y,a,weights=w).fit(cov_type='cluster',cov_kwds={'groups':block,'use_correction':True})
    # Independent least-squares and sandwich check, not a second statsmodels call.
    xx=a.T@(w[:,None]*a);xy=a.T@(w*y);beta=np.linalg.solve(xx,xy)
    inv=np.linalg.inv(xx);us=np.vstack([(a[block==g]*(w*(y-a@beta))[block==g,None]).sum(axis=0) for g in G])
    cov=inv@us.T@us@inv*ng/(ng-1)*(len(d)-1)/(len(d)-a.shape[1])
    check(f'{branch}_{ancestry}_{scope}_{model}_independent_beta',np.allclose(beta,result.params,rtol=1e-7,atol=1e-8))
    check(f'{branch}_{ancestry}_{scope}_{model}_independent_cov',np.allclose(cov,result.cov_params(),rtol=1e-6,atol=1e-8))
    se=math.sqrt(cov[1,1]);crit=t.ppf(.975,ng-1)
    row=dict(branch=branch,ancestry=ancestry,scope=scope,model=model,n=len(d),high_n=int(d.high.sum()) if 'high' in d else np.nan,blocks=ng,
             estimate=beta[1],cluster_se=se,ci_low=beta[1]-crit*se,ci_high=beta[1]+crit*se,
             p_pointwise=2*t.sf(abs(beta[1]/se),ng-1),condition_number=np.linalg.cond(a),
             bootstrap_valid=0,bootstrap_ci_low=np.nan,bootstrap_ci_high=np.nan,
             interpretation='post_result_conditional_association_not_selection_adjusted')
    if bootstrap or lobo:
        ax=np.stack([a[block==g].T@(w[block==g,None]*a[block==g]) for g in G]);ay=np.stack([a[block==g].T@(w[block==g]*y[block==g]) for g in G])
    if bootstrap:
        rng=np.random.default_rng(20260929);draw=[]
        for i in range(1000):
            counts=np.bincount(rng.integers(0,ng,ng),minlength=ng);X=np.einsum('g,gij->ij',counts,ax);Y=counts@ay
            if np.linalg.matrix_rank(X)==len(beta):draw.append(dict(iteration=i+1,estimate=np.linalg.solve(X,Y)[1],valid=True))
            else:draw.append(dict(iteration=i+1,estimate=np.nan,valid=False))
        v=pd.DataFrame(draw);row['bootstrap_valid']=int(v.valid.sum())
        if v.valid.sum()>=950:row['bootstrap_ci_low'],row['bootstrap_ci_high']=np.quantile(v.loc[v.valid,'estimate'],[.025,.975])
        write(v,f'bootstrap_{branch}_{ancestry}_{scope}_{model}.tsv')
    if lobo:
        rows=[]
        for i,g in enumerate(G):
            X=xx-ax[i];Y=xy-ay[i];valid=np.linalg.matrix_rank(X)==len(beta)
            b=np.linalg.solve(X,Y)[1] if valid else np.nan
            rows.append(dict(block=g,n_removed=int((block==g).sum()),high_removed=int(d.loc[block==g,'high'].sum()),estimate=b,delta=b-beta[1],valid=valid))
        write(pd.DataFrame(rows),f'lobo_{branch}_{ancestry}_{scope}_{model}.tsv')
    return row,result

def analyze():
    settings=json.loads((OUT/'liver_only_model_settings.json').read_text());rows=[];conting=[];sizes=[];distributions=[];matched=[]
    h={a:read(OUT/f'{a}_hcc_rank_scores.tsv').set_index('gene_symbol') for a in ['EUR','EAS']}
    for branch in ['boundary','neutral']:
        f=read(OUT/f'{branch}_liver_only_population.tsv')
        for ancestry in h:
            f[f'{ancestry}_Y']=f.gene_symbol.map(h[ancestry][f'{branch}_Y']);f[f'{ancestry}_fraction']=f.gene_symbol.map(h[ancestry][f'{branch}_rank_fraction'])
        write(f,f'{branch}_analysis_source.tsv')
        for ancestry in h:
            valid=f[f'{ancestry}_Y'].notna();nonseed=f.loc[~f.seed_gene&valid].copy()
            for threshold in [.01,.05,.1,.2]:
                nom=nonseed.liver_nomination_flag;aligned=nonseed[f'{ancestry}_fraction'].le(threshold)
                counts=[int((nom&aligned).sum()),int((nom&~aligned).sum()),int((~nom&aligned).sum()),int((~nom&~aligned).sum())]
                odds,p=fisher_exact(np.array(counts).reshape(2,2));conting.append(dict(branch=branch,ancestry=ancestry,threshold=threshold,nominee_aligned=counts[0],nominee_total=sum(counts[:2]),background_aligned=counts[2],background_total=sum(counts[2:]),OR=odds,p_descriptive=p))
            allpop=nonseed.loc[nonseed.complete_X&nonseed[f'{ancestry}_ld_block'].notna()].copy();allpop['Y']=allpop[f'{ancestry}_Y'];allpop['nominee']=allpop.liver_nomination_flag.astype(int)
            for name in ['overall_unadjusted','overall_X_adjusted']:
                x=pd.DataFrame({'intercept':1.,'nominee':allpop.nominee},index=allpop.index)
                if name.endswith('X_adjusted'):
                    for c in XCOL[1:]:x[c]=(allpop[c]-allpop[c].mean())/allpop[c].std(ddof=0)
                r,_=fit(allpop,x,ancestry,branch,'overall',name);rows.append(r)
            d=f.loc[f.support_population&valid&f[f'{ancestry}_ld_block'].notna()].copy();d['Y']=d[f'{ancestry}_Y']
            write(d,f'{branch}_{ancestry}_common_model_sample.tsv')
            sizes.append(dict(branch=branch,ancestry=ancestry,n=len(d),high=int(d.high.sum()),low=int(d.high.eq(0).sum()),high_blocks=d.loc[d.high.eq(1),f'{ancestry}_ld_block'].nunique(),blocks=d[f'{ancestry}_ld_block'].nunique()))
            for model in ['M0','M1','M2','M3']:
                x=design(d,model,settings[branch]);r,_=fit(d,x,ancestry,branch,'support',model,bootstrap=model in ['M2','M3'],lobo=model=='M3');rows.append(r)
            for scope,sub in [('S4',d[d.S_quintile.eq(4)]),('S5',d[d.S_quintile.eq(5)]),('exclude_B7plus',d[d.B.le(6)]),('EUR_EAS_common',d[d.EUR_Y.notna()&d.EAS_Y.notna()])]:
                if sub.high.nunique()<2:continue
                r,_=fit(sub,design(sub,'M3',settings[branch]),ancestry,branch,scope,'M3');rows.append(r)
            weights=1/d.groupby(f'{ancestry}_ld_block').gene_symbol.transform('size').to_numpy()
            r,_=fit(d,design(d,'M3',settings[branch]),ancestry,branch,'equal_block_weight','M3',weights=weights);rows.append(r)
            for b,g in d.groupby('B'):
                distributions.append(dict(branch=branch,ancestry=ancestry,B=b,n=len(g),blocks=g[f'{ancestry}_ld_block'].nunique(),raw_mean_Y=g.Y.mean(),raw_median_Y=g.Y.median(),q25=g.Y.quantile(.25),q75=g.Y.quantile(.75)))
            genes=read(BASE/f'27_locked_matching_design/{branch}_matched_liver_features.tsv').gene_symbol
            m=f.loc[f.gene_symbol.isin(genes)&valid&f[f'{ancestry}_ld_block'].notna()].copy();m['Y']=m[f'{ancestry}_Y']
            # Retain complete matched pairs after endpoint evaluability; no re-pairing.
            pairs=read(BASE/f'27_locked_matching_design/{branch}_matched_pairs.tsv')
            pairs=pairs[pairs.high_gene.isin(m.gene_symbol)&pairs.low_gene.isin(m.gene_symbol)]
            m=m[m.gene_symbol.isin(set(pairs.high_gene)|set(pairs.low_gene))]
            for model in ['M0','M2','M3']:
                x=design(m,model,settings[branch]);r,res=fit(m,x,ancestry,branch,'FAILED_match_diagnostic',model);rows.append(r)
                contribution=x.to_numpy()[:,2:]@res.params[2:] if model!='M0' else np.zeros(len(m))
                snp_columns=np.array([str(c).startswith('log_cis_snp_count') for c in x.columns])
                snp_contribution=x.to_numpy()[:,snp_columns]@res.params[snp_columns]
                matched.append(dict(branch=branch,ancestry=ancestry,model=model,pairs=len(pairs),n=len(m),estimate=r['estimate'],mean_fitted_covariate_difference=contribution[m.high.eq(1)].mean()-contribution[m.high.eq(0)].mean(),mean_fitted_SNP_component_difference=snp_contribution[m.high.eq(1)].mean()-snp_contribution[m.high.eq(0)].mean(),matching_gate='FAILED_not_reclassified'))
    write(pd.DataFrame(rows),'model_estimates.tsv');write(pd.DataFrame(conting),'overall_threshold_context.tsv');write(pd.DataFrame(sizes),'analysis_sample_counts.tsv');write(pd.DataFrame(distributions),'B_level_raw_distributions.tsv');write(pd.DataFrame(matched),'failed_match_model_diagnostic.tsv');write(pd.DataFrame(CHECKS),'model_validation.tsv')
    print(pd.DataFrame(rows).query('branch=="boundary" and scope=="support"').to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','scores','analyze']);args=p.parse_args();globals()[args.stage]()
