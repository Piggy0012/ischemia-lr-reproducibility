"""Generate the final ambient paragraph from audited numerical outputs."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parent;T=ROOT.parent/'outputs/reproducibility_v2/tables'

def read(name):return pd.read_csv(T/name,sep='\t',float_precision='round_trip')
def pct(x):return f'{100*float(x):.1f}%'

def main():
    model=json.loads((ROOT/'repro_decontx_model_narrative.json').read_text(encoding='utf-8'))
    assert model['status']=='complete'
    audit=json.loads((T/'ambient_rank_matched_audit.json').read_text(encoding='utf-8'))
    assert audit['status']=='complete' and audit['n_config_library_combinations']==22 and not audit['new_pvalues']
    assert audit['raw_and_corrected_fixed_universe']
    matched=read('ambient_rank_matched_summary.tsv');custom=read('ambient_concordance_summary.tsv')
    keys=read('ambient_rank_matched_candidate_keys.tsv')
    examples=read('ambient_fixed_examples.tsv')
    markers=read('ambient_markers_per_animal.tsv')
    markers=markers[(markers.config=='primary')&markers.gene.isin(['C1qa','Lyz2','Tyrobp'])].copy()
    assert len(markers)==66 and markers.raw_mean_count.gt(0).all()
    assert not markers.duplicated(['dataset','sample','cell_type','gene']).any()
    assert markers.corrected_mean_estimated_count.lt(markers.raw_mean_count).all()
    markers['reduction']=1-markers.corrected_mean_estimated_count/markers.raw_mean_count
    med=markers.groupby(['cell_type','gene']).reduction.median()
    assert len(med)==6 and markers.groupby(['cell_type','gene']).size().eq(11).all()
    def cr(config,gate,source='corrected',universe='common_corrected_eligible10'):
        d=custom[(custom.config==config)&(custom.comparison=='between_cohort_'+source)&
            (custom.universe==universe)&(custom.corrected_detection==gate)]
        assert len(d)==1;return d.iloc[0]
    def mr(config,source,metric):
        d=matched[(matched.config==config)&(matched.count_source==source)&(matched.metric==metric)&(matched.comparison=='between_cohorts')]
        assert len(d)==1;return d.iloc[0]
    for config in ['primary','reference_singlet']:
        k=keys[keys.config==config]
        assert len(k)>0 and not k.duplicated(['sender','receiver','ligand','receptor']).any()
        rows=matched[(matched.config==config)&(matched.comparison=='between_cohorts')]
        assert len(rows)==12 and rows.n_candidates.eq(len(k)).all()
        for row in rows.itertuples(index=False):
            assert abs(row.n_same_direction/row.n_candidates-row.same_direction_all_fraction)<1e-12
            if row.n_nonzero_both:
                assert abs(row.n_same_direction/row.n_nonzero_both-row.same_direction_nonzero_fraction)<1e-12
    # Keep the named examples auditable; fail instead of retaining stale prose.
    for ligand,receptor in [('Timp3','Kdr'),('Ptn','Ptprz1'),('Spp1','Itga5_Itgb1'),('Plat','Lrp1')]:
        e=examples[(examples.ligand==ligand)&(examples.receptor==receptor)]
        assert len(e)==4 and not e.duplicated(['dataset','config']).any()
        assert pd.api.types.is_bool_dtype(e.corrected_eligible10_ge1)
        assert (e.raw_score_difference*e.corrected_score_difference).gt(0).all()
        if ligand in ['Timp3','Ptn']:assert e.corrected_eligible10_ge1.all()
        elif ligand=='Spp1':assert not e.corrected_eligible10_ge1.any()
        else:
            lost=(e.dataset=='GSE174574')&(e.config=='primary')
            assert e.loc[lost,'corrected_eligible10_ge1'].eq(False).all()
            assert e.loc[~lost,'corrected_eligible10_ge1'].all()
    a,b=cr('primary','gt0'),cr('reference_singlet','gt0')
    c,d=cr('primary','ge1'),cr('reference_singlet','ge1')
    cr1=cr('primary','ge1','raw','common_raw_corrected_both_cohorts_eligible10')
    dr1=cr('reference_singlet','ge1','raw','common_raw_corrected_both_cohorts_eligible10')
    cc1=cr('primary','ge1','corrected','common_raw_corrected_both_cohorts_eligible10')
    dc1=cr('reference_singlet','ge1','corrected','common_raw_corrected_both_cohorts_eligible10')
    assert cr1.n_candidates==cc1.n_candidates==c.n_candidates
    assert dr1.n_candidates==dc1.n_candidates==d.n_candidates
    results_zh='### 校正后的髓系标记计数估计与候选资格敏感性\n\n'+model['results_zh']+'\n\n'
    results_zh+=f'原规则指定的目标细胞中，C1qa、Lyz2和Tyrobp的每文库每细胞平均估计计数均低于原始值；对六个细胞类型–基因组合分别汇总11个文库，文库降幅中位数的范围为{pct(med.min())}–{pct(med.max())}（图4A、B）。这些是模型扣减后的数值变化，不能独立验证其环境RNA来源。采用>0检测规则时，校正后共同候选为{int(a.n_candidates)}/{int(b.n_candidates)}条，原规则及参考支持的同向率为{int(a.n_same_direction)}/{int(a.n_candidates)}（{pct(a.same_direction_all_fraction)}，ρ={a.spearman_rho:.3f}）与{int(b.n_same_direction)}/{int(b.n_candidates)}（{pct(b.same_direction_all_fraction)}，ρ={b.spearman_rho:.3f}）。另设估计计数≥1的较严格检测阈值，并沿用相同的10%组别资格规则后，共同集合为{int(c.n_candidates)}/{int(d.n_candidates)}条，同向率为{int(c.n_same_direction)}/{int(c.n_candidates)}（{pct(c.same_direction_all_fraction)}）和{int(d.n_same_direction)}/{int(d.n_candidates)}（{pct(d.same_direction_all_fraction)}），ρ={c.spearman_rho:.3f}/{d.spearman_rho:.3f}（图4E、F）。在这两组较小的相同候选集合内，原始→校正的表达ρ分别为{cr1.spearman_rho:.3f}→{cc1.spearman_rho:.3f}和{dr1.spearman_rho:.3f}→{dc1.spearman_rho:.3f}。因此，较小集合的较高比例不能解释为对原全部候选的一致性改善。\n\n'
    results_en='### Corrected myeloid-marker estimates and eligibility sensitivity\n\n'+model['results_en']+'\n\n'
    results_en+=f'In primary-rule target cells, mean estimated counts per cell for C1qa, Lyz2 and Tyrobp were lower than raw values in every library. Each of the six cell-type–gene combinations was summarized across eleven libraries; median library-level reductions ranged from {pct(med.min())} to {pct(med.max())} (Figure 4A,B). These changes follow model-based subtraction and do not independently validate an ambient-RNA origin. Under the >0 detection rule, corrected common sets contained {int(a.n_candidates)}/{int(b.n_candidates)} candidates. Primary and reference-supported concordance was {int(a.n_same_direction)}/{int(a.n_candidates)} ({pct(a.same_direction_all_fraction)}; ρ={a.spearman_rho:.3f}) and {int(b.n_same_direction)}/{int(b.n_candidates)} ({pct(b.same_direction_all_fraction)}; ρ={b.spearman_rho:.3f}). Applying an additional, stricter estimated-count ≥1 detection threshold with the same 10% group eligibility rule yielded {int(c.n_candidates)}/{int(d.n_candidates)} common candidates, with concordance of {int(c.n_same_direction)}/{int(c.n_candidates)} ({pct(c.same_direction_all_fraction)}) and {int(d.n_same_direction)}/{int(d.n_candidates)} ({pct(d.same_direction_all_fraction)}) and ρ={c.spearman_rho:.3f}/{d.spearman_rho:.3f} (Figure 4E,F). Within these two smaller matched candidate sets, expression-effect correlations changed from {cr1.spearman_rho:.3f} to {cc1.spearman_rho:.3f} and from {dr1.spearman_rho:.3f} to {dc1.spearman_rho:.3f} from raw to corrected data. Higher agreement in the smaller sets therefore does not establish improvement across the original candidates.\n\n'
    p,q=mr('primary','decontx','custom_score'),mr('reference_singlet','decontx','custom_score')
    pn,qn=mr('primary','decontx','native_magnitude_priority'),mr('reference_singlet','decontx','native_magnitude_priority')
    pu,qu=mr('primary','decontx','diagnostic_unique_column_priority'),mr('reference_singlet','decontx','diagnostic_unique_column_priority')
    pr,qr=mr('primary','raw','custom_score'),mr('reference_singlet','raw','custom_score')
    pnr,qnr=mr('primary','raw','native_magnitude_priority'),mr('reference_singlet','raw','native_magnitude_priority')
    pur,qur=mr('primary','raw','diagnostic_unique_column_priority'),mr('reference_singlet','raw','diagnostic_unique_column_priority')
    results_zh+=f'为比较相同对象，将原始与校正数据、全部11个动物和全部六项度量均完整的目标集合固定为原规则{int(p.n_candidates)}条、参考支持{int(q.n_candidates)}条。原规则下，表达共可用性、原生优先级及独特分数列诊断的ρ（原始→校正）分别为{pr.spearman_rho:.3f}→{p.spearman_rho:.3f}、{pnr.spearman_rho:.3f}→{pn.spearman_rho:.3f}和{pur.spearman_rho:.3f}→{pu.spearman_rho:.3f}；参考支持下分别为{qr.spearman_rho:.3f}→{q.spearman_rho:.3f}、{qnr.spearman_rho:.3f}→{qn.spearman_rho:.3f}和{qur.spearman_rho:.3f}→{qu.spearman_rho:.3f}（图4C、D）。逐动物评分、全部及非零效应分母和全部度量结果同时存档，原始数据的置换P值不移用于这些描述性敏感性比较。Timp3–Kdr及Ptn–Ptprz1在两队列两种选择中保留≥1阈值下的资格，表达效应仍与原始数据同向；Spp1–Itga5/Itgb1在四种队列–选择组合均未达到该资格。Plat–Lrp1在原规则发现队列中未达到该资格，在其余三种组合中保留。未达到≥1资格不等于无表达或无通讯，也不能单凭模型确定这些转录本的细胞来源。'
    results_en+=f'To compare the same targets, candidates complete across raw/corrected data, all eleven animals and all six metrics were fixed at {int(p.n_candidates)} under primary selection and {int(q.n_candidates)} under reference-supported selection. Under primary selection, raw-to-corrected correlations were {pr.spearman_rho:.3f} to {p.spearman_rho:.3f} for coavailability, {pnr.spearman_rho:.3f} to {pn.spearman_rho:.3f} for native priority, and {pur.spearman_rho:.3f} to {pu.spearman_rho:.3f} for the unique-column diagnostic. Under reference-supported selection, the corresponding changes were {qr.spearman_rho:.3f} to {q.spearman_rho:.3f}, {qnr.spearman_rho:.3f} to {qn.spearman_rho:.3f}, and {qur.spearman_rho:.3f} to {qu.spearman_rho:.3f} (Figure 4C,D). Animal-level scores, all-candidate and nonzero-effect denominators, and every metric are retained; raw-data permutation P values do not apply to these descriptive sensitivity comparisons. Timp3–Kdr and Ptn–Ptprz1 retained eligibility under the ≥1 threshold in both cohorts and selections, with expression effects retaining their raw-data directions. Spp1–Itga5/Itgb1 did not meet that eligibility rule in any of the four cohort–selection combinations. Plat–Lrp1 did not meet it in primary-rule discovery cells but retained eligibility in the other three combinations. Failing ≥1 eligibility does not establish absence of expression or communication, nor can the model alone determine the cellular origin of these transcripts.'
    discussion_zh='环境RNA敏感性分析需要分别解释计数估计和候选资格。髓系标记计数的模型估计低于原始值，而固定候选集合中的相关变化因细胞选择与度量而异；这些计算结果不提供细胞来源的独立真值。单个原始UMI经估计扣减后可能低于1而仍为正，≥1阈值会改变低丰度条目的纳入范围，不能作为无表达的判据；>0资格保留也不能证明没有污染影响。DecontX以全脑参考标签分群，而下游保持既存细胞选择，参考不一致细胞中的真实表达可能被模型归为非本细胞来源组分。对Spp1等候选，应报告这一标签及检测阈值敏感性，而不据模型本身宣称其确定来源或通讯机制。'
    discussion_en='Ambient-RNA sensitivity requires separate interpretation of estimated counts and candidate eligibility. Myeloid-marker count estimates were lower than raw values, while correlations on fixed candidate sets changed differently across cell selections and metrics; these computations do not provide independent ground truth for cellular origin. A single observed UMI can become a positive estimate below one after subtraction, so the ≥1 threshold changes inclusion of low-abundance features and is not a criterion for absence of expression. Retained >0 eligibility likewise does not establish absence of contamination effects. DecontX groups cells by whole-brain reference labels while downstream cell selections remain fixed; genuine expression in reference-discordant cells may therefore be assigned to the non-native component. For candidates such as Spp1, label and detection-threshold sensitivity should be reported without claiming an established cellular source or communication mechanism from the model alone.'
    caption_zh='图4　环境RNA敏感性。A、B：原规则指定的星形胶质及内皮中C1qa、Lyz2、Tyrobp的每文库每细胞平均计数，纵轴作log1p转换；每条细线连接同一动物文库的原始值和DecontX估计值，共11个动物文库。C、D：原规则及参考支持选择，在raw/corrected、全部动物及全部六指标均完整的相同候选集合中比较跨队列效应ρ；蓝圆为原始值，橙三角为校正估计。Native和Unique分别为原生及独特分数列诊断优先级。E、F：两队列共同通过10%组别门槛的候选数量；原始>0、校正>0与校正≥1对应不同检测定义，集合可不同。≥1为另设的敏感性阈值，未达门槛不证明无表达。图中没有新增显著性检验；计数估计、检测阈值与数值收敛来源选择均单独记录。'
    caption_en='Figure 4. Ambient-RNA sensitivity. A,B: Log1p-transformed library-level mean counts per cell for C1qa, Lyz2 and Tyrobp in primary-rule astrocytes and endothelium. Thin lines connect raw and DecontX-estimated values from the same animal library (11 libraries). C,D: Primary and reference-supported selections compare cross-cohort effect correlations on the same candidates complete for raw/corrected data, every animal and all six metrics. Blue circles denote raw values and orange triangles corrected estimates. Native and Unique denote native and unique-column diagnostic priorities. E,F: Common candidate counts passing the 10% group gate under raw >0, corrected >0 and corrected ≥1 detection definitions; these eligible sets may differ. The ≥1 rule is an additional sensitivity threshold, and failure does not establish absence of expression. No new inferential tests are shown. Estimated counts, detection thresholds and numerical-convergence source selection are documented separately.'
    out={'status':'complete','methods_zh':model['methods_zh'],'methods_en':model['methods_en'],
        'results_zh':results_zh,'results_en':results_en,'discussion_zh':discussion_zh,'discussion_en':discussion_en,
        'caption_zh':caption_zh,'caption_en':caption_en}
    (ROOT/'repro_ambient_narrative.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Ambient narrative generated from complete matched results')

if __name__=='__main__':main()
