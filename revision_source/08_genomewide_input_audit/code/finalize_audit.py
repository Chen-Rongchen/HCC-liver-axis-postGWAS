"""Close the bounded input audit with unchanged-source checks and a readout."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from audit_inputs import BASE,ROOT,hashes

RES=BASE/'results';checks=[];cache={}
def sha(p):
    if p not in cache:cache[p]=hashes(p)
    return cache[p]['sha256']
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
def check(name,ok,detail=None):
    checks.append(dict(name=name,passed=bool(ok),detail=detail))
    if not ok:raise AssertionError(checks[-1])

collision=json.loads((RES/'collision_tests.json').read_text())
replay=json.loads((RES/'score_replay_checks.json').read_text())
check('actual_legacy_collision_regressions',collision['all_passed'] and all(x['passed'] for x in collision['checks']))
check('full_gene_replay_checks',replay['all_passed'] and all(x['passed'] for x in replay['checks']))
model=read(RES/'M2_impact.tsv');enrich=read(RES/'enrichment_impact.tsv');members=read(RES/'membership_impact.tsv')
archived=read(ROOT/'analysis_hbsn_restructure/35_seed_region_masking/out_06_rule_A_outcomes/M2.tsv')
for r in model[model.version.eq('frozen')].itertuples():
    old=archived[archived.ancestry.eq(r.ancestry)&archived.part.eq(r.part)].iloc[0]
    check(f'{r.ancestry}_{r.part}_existing_M2_replay',np.allclose([r.estimate,r.cluster_se,r.ci_low,r.ci_high],[old.estimate,old.cluster_se,old.ci_low,old.ci_high],rtol=1e-10,atol=1e-14))
for (anc,part),g in enrich.groupby(['ancestry','part']):
    check(f'{anc}_{part}_enrichment_counts_preserved',all(g[c].nunique()==1 for c in ['nominee_top5','nominee_total','background_top5','background_total']))
check('candidate_members_preserved',members.candidate_added.eq(0).all() and members.candidate_removed.eq(0).all())
summaries={a:json.loads((RES/a/'summary.json').read_text()) for a in ['GCST90860790','GCST90860791']}
all_sources=pd.concat([read(RES/a/'source_hashes.tsv') for a in summaries]+[read(RES/'score_replay_sources.tsv')]).drop_duplicates('path')
check('all_raw_cache_and_mapping_inputs_unchanged',all(sha(ROOT/r.path)==r.sha256 for r in all_sources.itertuples()),len(all_sources))
protected=json.loads((ROOT/'analysis_signal_resolution_c/00_registry/pre_execution_hashes.json').read_text())
check('all_138_protected_files_unchanged',len(protected)==138 and all(sha(ROOT/p)==h for p,h in protected.items()),len(protected))
check('protocol_unchanged',sha(BASE/'PROTOCOL.txt')==json.loads((BASE/'protocol_lock.json').read_text())['sha256'])
check('no_C2_models',not any((ROOT/'analysis_signal_resolution_c/04_matched_comparison').iterdir()))
events=read(ROOT/'github_release/HCC-liver-axis-postGWAS/derived/support_events.tsv')
for a in summaries:
    genes=read(RES/a/'P_cache_affected_genes.tsv')
    check(a+'_conflict_affected_genes_outside_283_candidates',not set(genes.gene)&set(events.gene_symbol))

delta=model.pivot(index=['ancestry','part'],columns='version',values='estimate')
report=dict(status='bounded_input_audit_complete_deterministic_collision_repaired',
    raw_rows={a:s['counters']['raw_rows'] for a,s in summaries.items()},
    conflicting_keys={'EUR':48,'EAS':1},boundary_scores_affected={'EUR':61,'EAS':61},
    conservative_roundtrip_exclusion_keys={'EUR':1322,'EAS':970},
    mapping_note='Forward unique but reverse-chain ambiguous; exclusion sensitivity does not establish physical variant identity.',
    new_scientific_model_specifications=0,existing_M2_fits=replay['existing_M2_refits'],
    max_absolute_collision_repair_M2_change=float((delta.collision_repaired-delta.frozen).abs().max()),
    max_absolute_roundtrip_sensitivity_M2_change=float((delta.roundtrip_exclusion_sensitivity-delta.frozen).abs().max()),
    C1_formal_admissions=0,C4='Stop—technical',new_shared_signal_evidence=False,
    original_cache_overwritten=False,manuscript_updated=False,submission_updated=False,
    submission_sha256=sha(ROOT/'submission/EJHG_Review_Package.zip'),
    collision_checks=len(collision['checks']),score_replay_checks=len(replay['checks']),closeout_checks=checks)
(RES/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
text='''2026-10-04 全基因组输入核查与有限修复

结论
已完成 EUR/EAS 原始 HCC 输入到全部冻结缓存的核查，并发现、定位和修复一个真实适配器缺陷。
不同 GRCh38 位点映射到同一 GRCh37 无序等位基因键，若出现在同一个 chunk，旧适配器会把二者都视为 fresh，NumPy 重复索引赋值令后值覆盖前值，漏记冲突。
涉及 EUR 48 键（chr1:7、chr6:2、chr13:2、chr14:16、chr15:21）、EAS 1 键（chr19）。现有 7 个缓存字段合计差异为 336 和 7，不是 343 个独立变异。
逐个使用原始适配器方法在内存中重放全部 49 个冲突键，重现了覆盖结果；修复合并器对行顺序、同批/跨批、相同重复、冲突后追加均有回归检查。

修复方式
不改写原始 GWAS 或原缓存。每个终点提供 collision_repair_overlay.tsv；code/collision_repair.py 提供带原缓存哈希检查的加载器和永久记录冲突的合并器。修复后所有冲突键的 p/beta/se/eaf/n 置缺失、effect_token=-1、conflict=1。
完整扫描器 audit_inputs.py 也已在同批和跨批级别处理冲突，恢复与既定冲突排除契约一致的结果。旧 qc_primary_2026_endpoints.py 保留作历史，不应作为新稿的生成入口。

扫描范围
EUR 21,986,982 行，EAS 7,735,435 行；原始文件 MD5 与冻结的官方校验值一致。
每个终点检查 22 染色体 × 7 数值/状态字段。EUR 301,087 行不满足历史数值契约，但这些行没有精确匹配冻结的坐标/等位基因键；EAS 无此类行。不能据此推广为所有来源数值语义已经核验。
在保留的同染色体正链转换范围内，未发现正向多重命中、唯一正向映射遗漏或坐标差异。
另有 EUR 1,370、EAS 971 条原始匹配记录处于回映射歧义区间，对应 1,322、970 个键，其中包括上述冲突键。这与正向一对多不是同一概念。仅凭 chain 不为它们猜测原始物理身份。

基因分数、排名与成员
对 18,321 个冻结基因，独立重放 boundary/neutral 稳定 ACAT、P=1、排名分数和 Y，与归档一致；保留重复的映射记录和原排名总体。16 个合成边界用例另与 400 位精度结果对照。
冲突修复使 EUR/EAS 各 61 个基因的分数发生变化，这些基因均不在当前 283 个候选内。Top-5% 集合、283 个 EUR 候选、1,652 个共同模型成员及规则 A 后 1,259 个成员不变。
严格排除所有回映射歧义键另列为 roundtrip_exclusion_sensitivity：EUR/EAS 各 805/756 个基因分数变化，但上述 top-5%、候选和模型成员仍不变。这是保守敏感性，不冒充来源身份修复。

既有结果的必要重算（无新增模型规格）
                       原冻结值       冲突修复后        回映射排除敏感性
EUR 原始 M2            0.691839        0.691846          0.691892
EUR 规则 A 后 M2       0.002534        0.002541          0.002655
EAS 原始 M2            0.255712        0.255711          0.255644
EAS 规则 A 后 M2       0.155119        0.155116          0.154922
完整 SE、区间、分母见 results/M2_impact.tsv。仅冲突修复最大点估计变化约 6.95e-6；更保守排除最大变化约 1.96e-4。
EUR 描述性富集 OR 2.6274→1.5479、EAS 2.0514→1.7442，所有列联表计数和 Fisher P 不变。
283 候选未变，且肝脏性状输入未修改，本次没有重算支持事件配对或共定位，不把这一点写成独立重建 544 配对。

科学意义与边界
这次是确定性代码/输入修复及其影响核查，不是获得新的疾病遗传学补强。未改变“EUR 广度条件估计对规则 A 的成员组成敏感”这一判断，也没有证明 input reuse 因果性地解释 HCC 关联。
规则 A 改变估计总体；0.692→0.003 并非在同一总体中识别出某种介导效应，后者区间仍宽，不能称为证明广度没有作用。
HCC 有符号效应、Beta 尺度、N 语义以及 Track C 的 LD/覆盖问题仍按已有来源证据状态处理。C1 仍 0/19，C4 Stop—technical；没有运行 C2、Genomic SEM 或新增共定位。
全基因组核查仅覆盖 HCC EUR/EAS 原始输入、冻结映射缓存以及现有分数/排名/模型链；没有独立复现全部 11 个肝脏性状原始来源。

投稿定位
可以据此讨论以 cross-trait evidence interpretation / regional genetic representation 重写 Track A。相同 gene 名字不等于相同科学结果，但新的表述仍需与文献逐项比较，不能仅靠改标题证明新颖性。
规则 A 应进入主结果，source cis 聚合排名不等于显著 eQTL/效应基因证据，known anchors 不称新发现基因。上述修复及来源语义限制也需要进入新稿方法与公开复现入口。
尚未改稿、改项目计划或投稿包；真实作者声明仍待作者提供。唯一 EJHG ZIP 与 138 个受保护文件保持原哈希。

复核入口
results/validation.json：最终状态与保护检查。
results/score_replay_checks.json：全基因重放。
results/collision_tests.json：真实旧代码缺陷重现及修复回归。
results/GCST*/conflicting_source_records.tsv：逐条原始冲突记录。
results/GCST*/mapping_and_numeric_exceptions.tsv.gz：完整回映射异常记录。
results/GCST*/source_hashes.tsv 与 results/score_replay_sources.tsv：输入版本。
AUDIT_MANIFEST.tsv：本次脚本和输出清单。
运行顺序：audit_inputs.py GCST90860790、audit_inputs.py GCST90860791；replay_scores.py；test_collisions.py；finalize_audit.py。
扫描器拒绝复用已有扫描结果目录；再次完整扫描应使用独立审计目录，以免异常日志重复。该重跑保护在首次扫描完成后加入，不改变已生成数据。不得直接将旧缓存或旧稿作为修复后的最终版。
'''
(BASE/'AUDIT_REPORT.txt').write_text(text)
entry=ROOT/'analysis_signal_resolution_c/README.md';old=entry.read_text()
marker='2026-10-04 全基因组输入核查已完成。'
if marker not in old:
    old=old.replace('# Track C 当前入口与历史记录\n','# Track C 当前入口与历史记录\n\n'+marker+'见 [核查与有限修复](08_genomewide_input_audit/AUDIT_REPORT.txt) 和 [结果](08_genomewide_input_audit/results/validation.json)。新增发现 EUR 48/EAS 1 个同批映射冲突漏检；已提供只读缓存修复层并重放既有结果，三位小数结论和候选成员不变。该更新不解除 C1 技术限制；原包未改。\n',1)
    entry.write_text(old)
manifest=[]
for p in sorted(BASE.rglob('*')):
    if p.is_file() and '__pycache__' not in p.parts and p.name!='AUDIT_MANIFEST.tsv':manifest.append({**hashes(p),'path':str(p.relative_to(BASE))})
pd.DataFrame(manifest).to_csv(BASE/'AUDIT_MANIFEST.tsv',sep='\t',index=False)
branch=ROOT/'analysis_signal_resolution_c';inventory=[]
for p in sorted(branch.rglob('*')):
    if p.is_file() and '__pycache__' not in p.parts and p.name!='ARTIFACT_MANIFEST.tsv' and p.suffix!='.log':inventory.append({**hashes(p),'path':str(p.relative_to(branch))})
pd.DataFrame(inventory).to_csv(branch/'ARTIFACT_MANIFEST.tsv',sep='\t',index=False)
print(json.dumps({k:v for k,v in report.items() if k!='closeout_checks'},indent=2,ensure_ascii=False))
