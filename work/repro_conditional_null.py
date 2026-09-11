"""Explain joint versus single-cohort label-null evidence from saved allocations."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2/tables'
STATS=['same_direction_all_fraction','same_direction_nonzero_fraction','spearman_rho']
TOL=1e-12

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    plan=ROOT/'repro_conditional_null_plan.md'
    sources={'fixed_common':('null_fixed_common_all_allocations.tsv.gz','null_fixed_common_summary.tsv','A_fixed_common_metrics'),
             'dynamic_gate':('null_dynamic_gate_all_allocations.tsv.gz','null_dynamic_gate_summary.tsv','C_exploratory_dynamic_gate')}
    all_values=[];all_summaries=[];hashes={str(plan):sha(plan)};checks=[]
    for analysis,(name,summary_name,family) in sources.items():
        path=OUT/name;sp=OUT/summary_name
        hashes[name]=sha(path);hashes[summary_name]=sha(sp)
        data=pd.read_csv(path,sep='\t');joint=pd.read_csv(sp,sep='\t')
        for config,block in data.groupby('config',sort=False):
            obs=block[block.is_observed]
            od=obs.discovery_allocation.unique();ov=obs.validation_allocation.unique()
            assert len(od)==len(ov)==1
            for randomized,fixed,filter_column,observed_value,n in [
                ('GSE174574','GSE245386','validation_allocation',ov[0],20),
                ('GSE245386','GSE174574','discovery_allocation',od[0],10)]:
                subset=block[block[filter_column]==observed_value].copy()
                subset.insert(0,'fixed_observed_cohort',fixed)
                subset.insert(0,'randomized_cohort',randomized)
                subset.insert(0,'analysis',analysis)
                subset['n_conditional_allocations']=n
                all_values.append(subset)
                for metric,b in subset.groupby('metric',sort=False):
                    assert len(b)==n and b.is_observed.sum()==1
                    for stat in STATS:
                        observed=float(b.loc[b.is_observed,stat].iloc[0])
                        values=b[stat].to_numpy(float)
                        assert np.isfinite(values).all()
                        k=int((values>=observed-TOL).sum());p=k/n
                        j=joint[(joint.config==config)&(joint.metric==metric)&(joint.statistic==stat)]
                        assert len(j)==1 and abs(j.observed.iloc[0]-observed)<=TOL
                        assert p>=1/n and np.isclose(p*n,round(p*n))
                        all_summaries.append({'analysis':analysis,'config':config,'metric':metric,'statistic':stat,
                            'randomized_cohort':randomized,'fixed_observed_cohort':fixed,
                            'observed':observed,'null_mean':float(values.mean()),'null_median':float(np.median(values)),
                            'null_q025':float(np.quantile(values,.025)),'null_q975':float(np.quantile(values,.975)),
                            'n_conditional_allocations':n,'n_at_least_observed':k,
                            'conditional_unadjusted_upper_tail_p':p,'minimum_attainable_upper_tail_p':1/n,
                            'joint_n_allocations':200,'joint_exact_upper_tail_p':float(j.exact_upper_tail_p.iloc[0]),
                            'joint_holm_adjusted_p':float(j.holm_adjusted_p.iloc[0]),'original_joint_test_family':family,
                            'observed_n_candidates':int(b.loc[b.is_observed,'n_candidates'].iloc[0]),
                            'conditional_n_candidates_min':int(b.n_candidates.min()),
                            'conditional_n_candidates_max':int(b.n_candidates.max()),
                            'conditional_p_adjusted_for_multiple_testing':False})
                        checks.append({'analysis':analysis,'config':config,'metric':metric,'statistic':stat,
                            'randomized_cohort':randomized,'passed':True})
    values=pd.concat(all_values,ignore_index=True);summary=pd.DataFrame(all_summaries)
    assert len(summary)==72
    values.to_csv(OUT/'null_conditional_all_allocations.tsv.gz',sep='\t',index=False)
    summary.to_csv(OUT/'null_conditional_summary.tsv',sep='\t',index=False)
    plan_text=plan.read_text(encoding='utf-8')
    (OUT/'null_conditional_methods.md').write_text(plan_text,encoding='utf-8')
    explanation='在已知联合零分布结果后，追加条件性解释诊断：分别固定另一队列的真实处理标签，仅穷举发现队列的20种或外部队列的10种标签分配，完整复用既存多变量候选分数与逐配置表达门槛。单侧条件精确P为对应配置中统计量不小于观测值的比例，保留原配置；这些P仅作解释性诊断，不修改原A、B、C校正家族，也不构成新增确认性检验。两种细胞选择的动态10%门槛同向率及ρ在联合200配置中均得到P=0.005，但只重排发现队列时为0.05，只重排外部队列时为0.10。联合随机化评估两研究标签均无关联的整体零设定；拒绝它不足以证明两个队列分别均具有可重复的疾病关联。若要建立这种复制主张，需要另行论证针对“至少一个队列无相关关联”的联合原假设的检验。本诊断不将两个条件P的最大值包装为已验证的新检验；外部仅2个缺血动物产生的分辨率与设计限制仍然存在。'
    # The numeric statement is verified against every dynamic statistic/config.
    dyn=summary[summary.analysis=='dynamic_gate']
    assert np.allclose(dyn[dyn.randomized_cohort=='GSE174574'].conditional_unadjusted_upper_tail_p,.05)
    assert np.allclose(dyn[dyn.randomized_cohort=='GSE245386'].conditional_unadjusted_upper_tail_p,.10)
    (OUT/'null_conditional_explanation_zh.md').write_text(explanation+'\n',encoding='utf-8')
    audit={'status':'complete','new_confirmatory_test':False,'original_families_changed':False,
        'joint_results_known_before_diagnostic_plan':True,'source_sha256':hashes,'script_sha256':sha(Path(__file__)),
        'n_checks':len(checks),'checks':checks,'no_maximum_pvalue_replication_test_claim':True,
        'conditional_pvalues_are_unadjusted_diagnostics':True,'distribution_rows':len(values),'summary_rows':len(summary)}
    (OUT/'null_conditional_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    display=summary[(summary.metric.isin(['custom_score','magnitude_priority']))&
                    (summary.statistic.isin(['same_direction_all_fraction','spearman_rho']))]
    print(display[['analysis','config','metric','statistic','randomized_cohort','observed',
                   'conditional_unadjusted_upper_tail_p','joint_exact_upper_tail_p']].to_string(index=False))

if __name__=='__main__':main()
