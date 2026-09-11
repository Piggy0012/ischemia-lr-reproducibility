"""Prepare source-bound bilingual CellChat extension prose after independent audit."""
from pathlib import Path
from datetime import datetime, timezone
import json
import hashlib

ROOT = Path(__file__).resolve().parent
audit_path = ROOT / 'repro_v3_cellchat_editorial_audit.json'
review = json.loads(audit_path.read_text(encoding='utf-8'))
assert review['status'] == 'complete'
assert review['all_11_zero_strength_targets'] == 3
assert review['all_11_zero_strength_targets_with_nonzero_rank_effect'] == 3
assert [x['n_same_direction'] for x in review['primary_metrics']] == [23, 10, 16, 13]

narrative = {
    'status': 'complete',
    'scope': 'exploratory controlled-resource native CellChat extension; descriptive results only',
    'created_at_utc': datetime.now(timezone.utc).isoformat(),
    'main_methods_en': (
        'We added a controlled-resource analysis with the native CellChat R package '
        '(2.2.0.9001) using the same 11 raw-count libraries and frozen primary cell labels. '
        'The catalogue contained 1,548 shared, gene-covered interaction definitions with native '
        'complexes and cofactors. Rankings were calculated on 1,222 source–target–interaction '
        'entries with complete archived LIANA scores in every library. We compared native '
        'CellChat inferred strength, an analyst-derived strength percentile, expression '
        'coavailability and LIANA RRA recomputed once per distinct score column on this fixed '
        'network. These additional cohort contrasts were descriptive; the original permutation '
        'P values were not applied to them. [@Jin2025CellChat][@Jin2021CellChat]'),
    'main_methods_zh': (
        '使用原生 CellChat R 包（2.2.0.9001），对相同11个原始计数文库及固定的主要细胞标签'
        '进行了控制资源范围的扩展分析。目录包含1,548个共享且基因覆盖完整的互作定义，保留'
        '原生复合物和辅助因子。排名背景固定为在所有文库中均具有完整既存 LIANA 得分的'
        '1,222个来源细胞—靶细胞—互作条目。比较指标包括 CellChat 原生推断强度、研究者'
        '计算的强度百分位、表达共可用性，以及在该固定网络中每个不同得分列仅排名一次后'
        '重算的 LIANA RRA。新增队列对比仅作描述性分析，未套用原置换检验的P值。'
        '[@Jin2025CellChat][@Jin2021CellChat]'),
    'main_results_en': (
        'Exact native interaction definitions retained 32 of the original 187 directed '
        'candidates (29 distinct ligand–receptor definitions; 8 astrocyte→endothelial and 24 '
        'reverse-direction entries). On these same 32 entries, coavailability was concordant '
        'for 23/32 (ρ=0.565), native CellChat strength for 10/32 (ρ=0.549), the derived strength '
        'percentile for 16/32 (ρ=0.364), and matched-network LIANA RRA priority for 13/32 '
        '(ρ=0.451). Native strength had zero cohort effects for 19 entries in discovery and '
        '3 externally; concordance among the 13 entries with nonzero effects in both cohorts '
        'was 10/13. The corresponding nonzero denominator for LIANA RRA was 19 (13/19 '
        'concordant). Three entries had zero native strength in all 11 libraries but nonzero '
        'percentile contrasts because the zero-strength tie block changed across libraries '
        '(Supplementary Figure S6).'),
    'main_results_zh': (
        '原生互作定义的精确匹配保留了原187个有向候选中的32个条目，涉及29个不同的'
        '配体—受体定义，其中星形胶质细胞→内皮细胞8条、反向24条。在相同32条中，'
        '表达共可用性同向23/32（ρ=0.565），CellChat 原生强度同向10/32（ρ=0.549），'
        '研究者计算的强度百分位同向16/32（ρ=0.364），相同网络中重算的 LIANA RRA '
        '优先级同向13/32（ρ=0.451）。原生强度的队列效应在发现队列有19条为零，'
        '外部队列有3条为零；两队列效应均非零的13条中，同向10/13。LIANA RRA '
        '对应的非零分母为19条，同向13/19。另有3条在全部11库中的原生强度均为零，'
        '仍因不同文库中零强度并列块的大小变化而具有非零百分位对比（补充图S6）。'),
    'main_discussion_en': (
        'Native inferred strength and its derived percentile describe different features '
        'of a candidate. In this controlled extension, zero strengths reduced the number '
        'of nonzero cohort contrasts, while the percentile position of the zero-strength '
        'tie block shifted with the network background. Together with the resource '
        'attrition, these results support preserving native strengths, the exact ranking '
        'background and zero-block diagnostics when comparing candidate priorities '
        'across cohorts within an explicit resource and cell-selection scope.'),
    'main_discussion_zh': (
        '原生推断强度及其衍生百分位描述候选的不同特征。在本次控制范围的扩展中，'
        '零强度缩小了非零队列效应的比较集合，而零强度并列块的百分位位置随网络'
        '背景发生移动。结合资源匹配导致的候选减少，这些结果支持在明确资源和'
        '细胞选择范围的前提下，同步保留原生强度、确切的排名背景及零并列块诊断，'
        '以便比较候选优先级的跨队列变化。'),
    'supplement_methods_en': [
        'The exploratory extension used the complete native CellChat R package 2.2.0.9001 '
        '(official upstream commit 75253cd0c9e68410e6e721a6d3a0419a1d7e358f) in R 4.6.1. '
        'The Windows binary was a third-party r-universe build of the official source; its '
        'build provenance, installed core functions and complete native mouse database were '
        'verified against the pinned source. Each of the original 11 libraries was analyzed '
        'separately: GSE174574 contained three Sham and three MCAO libraries, and GSE245386 '
        'contained three Sham and two MCAO libraries. Frozen primary QC labels and all '
        'context cell types with at least 30 cells were retained. Third-cohort, reference-selected '
        'and DecontX-corrected inputs were outside this extension. [@Jin2025CellChat][@Jin2021CellChat]',
        'Complete complex subunits were expanded and canonically ordered without orthology '
        'conversion or alias substitution. Among native mouse protein-signaling pairs shared '
        'with mouseconsensus, 1,548 had unambiguous native rows and complete ligand, receptor '
        'and referenced cofactor gene coverage in all original matrices. Native complex and '
        'cofactor definitions were retained. Per-cell log1p(10,000×count/all-gene UMI total) '
        'values reproduced the archived LIANA float32 normalization, then were losslessly '
        'promoted to float64 for R input. Signaling genes were exported only after normalization; '
        'their subset totals were not used as denominators.',
        'Each native CellChat object used this controlled catalogue, subsetData and '
        'computeCommunProb with explicit LR.use, type="triMean", raw.use=TRUE, '
        'population.size=FALSE, nboot=100, seed.use=20260911, Kh=0.5 and n=1. '
        'Execution was sequential, with no PPI projection or spatial-distance modeling. '
        'Explicit LR.use bypassed per-library overexpressed-interaction preselection; thus '
        'this was a controlled configuration rather than the default end-to-end workflow. '
        'The native routine additionally scales signaling expression by its library-specific '
        'maximum. Its net$prob is termed probability by the package but represents inferred '
        'communication strength, not a calibrated biological probability. All scores, including '
        'zeros, and internal cell-label P values were retained. These internal P values neither '
        'filtered this comparison nor served as animal-level disease tests.',
        'The comparison network was fixed before reading new cross-cohort CellChat results: '
        '1,222 source–target–interaction keys had complete LIANA component scores in all '
        '11 libraries on the controlled resource and were present in every native output. '
        'This selection conditions the comparison on LIANA detection eligibility. For each '
        'library, the analyst-derived CellChat percentile was (ascending average rank of '
        'native strength−1)/(1,222−1), so larger values indicate higher relative strength. '
        'This is not a native CellChat consensus rank. LIANA unique-column RRA was '
        'recomputed over precisely the same 1,222 entries using lr_means, expr_prod and '
        'lrscore once each, with the archived float32 component representation restored. '
        'The RRA priority was 1−RRA. Original ranks from other network backgrounds were '
        'not reused. The original expression coavailability formula was evaluated on the '
        'same 32 target entries.',
        'For each metric, the cohort contrast was the MCAO mean minus the Sham mean across '
        'library-level scores. Spearman correlation, same-direction counts over all 32 entries, '
        'same-direction counts among entries nonzero in both cohorts, and zero-effect counts '
        'were reported separately. The existing absolute effect tolerance of 10^−12 was applied '
        'and audited. Native RDS scores were preserved; lossless 17-significant-digit exports '
        'were verified against them before round-trip parsing. Original A/B/C permutation '
        'P values do not apply to the changed metrics and candidate set, and no new '
        'inferential disease tests were run. Direction and cofactor subsets were descriptive '
        'only. The separate full native-resource sensitivity arm was not run; this analysis '
        'does not assess CellChat performance across its complete native catalogue.'
    ],
    'supplement_methods_zh': [
        '该探索性扩展在R 4.6.1中使用完整原生 CellChat R 包2.2.0.9001，固定官方上游提交为'
        '75253cd0c9e68410e6e721a6d3a0419a1d7e358f。Windows二进制包来自 r-universe '
        '对官方源码的第三方构建；构建来源、已安装核心函数及完整原生小鼠数据库均与固定'
        '源码进行了核对。原11库分别分析：GSE174574为3个Sham与3个MCAO文库，'
        'GSE245386为3个Sham与2个MCAO文库。采用固定的主要QC细胞标签，并保留每库'
        '至少30个细胞的全部背景细胞类型。本扩展未纳入第三队列、参考筛选或DecontX'
        '校正后的输入。[@Jin2025CellChat][@Jin2021CellChat]',
        '展开完整复合物亚基并按规范顺序排列，未作跨物种转换或别名替换。在原生小鼠'
        '蛋白信号目录与mouseconsensus共享的配体—受体定义中，1,548项具有唯一的原生行'
        '映射，且其配体、受体和所引用辅助因子基因在全部原始矩阵中均有覆盖。保留原生'
        '复合物和辅助因子定义。每细胞按log1p(10,000×计数/全部基因UMI总量)标准化，'
        '复现既存LIANA的float32运算后，无损提升至float64作为R输入。仅在完成标准化'
        '后导出信号基因，未使用该基因子集的计数总量作为分母。',
        '每库创建原生CellChat对象并指定控制目录，调用subsetData后，以显式LR.use、'
        'type="triMean"、raw.use=TRUE、population.size=FALSE、nboot=100、'
        'seed.use=20260911、Kh=0.5和n=1调用computeCommunProb。逐库顺序执行，'
        '不作PPI投影或空间距离建模。显式LR.use绕过了每库过表达互作预筛选，因此属于'
        '控制配置，非默认的完整分析流程。原生函数还按每库信号表达最大值进行缩放。'
        '包中的net$prob虽命名为probability，实为推断通讯强度，并非校准后的生物学'
        '概率。保留全部零得分及库内细胞标签检验P值；这些P值既不用于筛选本次比较，'
        '也不作为动物层面的疾病检验。',
        '在读取新增跨队列CellChat结果前，比较网络固定为控制资源中、在所有11库均有'
        '完整LIANA分量得分的1,222个来源细胞—靶细胞—互作条目；这些条目也均存在于'
        '各原生输出。该选择使比较以LIANA检测资格为条件。每库由研究者计算的强度'
        '百分位为（原生强度的升序平均秩−1）/(1,222−1)，数值越大表示相对强度越高，'
        '并非CellChat原生共识排名。在完全相同的1,222条上重算LIANA单列RRA，'
        'lr_means、expr_prod及lrscore各使用一次，并恢复既存分量的float32数值表示。'
        'RRA优先级定义为1−RRA，未复用其他网络背景的既存排名。相同32个目标条目'
        '也按原公式计算表达共可用性。',
        '各指标的队列对比均为文库层面MCAO均值减Sham均值。分别报告Spearman相关系数、'
        '以全部32条为分母的同向计数、以两队列效应均非零条目为分母的同向计数，以及'
        '零效应数。应用并审计原有10^−12绝对效应容差。保留原生RDS，通过17位有效数字'
        '的无损导出与原RDS核对后，再以可往返精度读取。原A/B/C置换检验P值不适用于'
        '改变后的指标和候选集合，本扩展未进行新的疾病效应推断检验。方向与辅助因子'
        '分层仅作描述。单独的完整原生资源敏感性分析未运行，因此本分析不评估CellChat'
        '在其完整原生目录中的性能。'
    ],
    'supplement_results_en': [
        'The resource audit retained 32/187 original directed candidates (17.1%): 8 '
        'astrocyte→endothelial and 24 endothelial→astrocyte entries, representing 29 '
        'distinct ligand–receptor definitions. The remaining 155 lacked an exact complete '
        'subunit-pair match in native CellChatDB; none was additionally lost at the common '
        'full-network intersection. These exclusions concern resource definitions, not '
        'failed CellChat predictions. Thirty-one retained entries lacked native cofactors '
        'and one had a cofactor, precluding a meaningful cofactor-stratum correlation.',
        'On the common 32 entries, coavailability concordance was 23/32 (71.9%; ρ=0.565); '
        'native strength concordance was 10/32 (31.3%; ρ=0.549); derived percentile '
        'concordance was 16/32 (50.0%; ρ=0.364); and matched-network LIANA RRA priority '
        'concordance was 13/32 (40.6%; ρ=0.451). Native strength had 19 zero effects in '
        'discovery and 3 externally, leaving 13 entries nonzero in both cohorts; 10/13 '
        '(76.9%) were concordant. LIANA RRA had 10 and 12 zero effects, leaving 19 '
        'nonzero-both entries; 13/19 (68.4%) were concordant. Coavailability and the '
        'derived percentile had no zero cohort effects. The 10^−12 tolerance removed '
        'no nonzero contrasts for these four metrics in either cohort.',
        'Native strength was zero in all 11 libraries for astrocyte→endothelial '
        'Lamb2–Dag1 and Lamb2–Itga6_Itgb1, and endothelial→astrocyte Entpd1–Adora2b. '
        'Nevertheless, each acquired a derived percentile contrast of approximately '
        '+0.014060 in discovery and −0.024775 externally. Across libraries, 160–862 '
        'of the 1,222 network strengths were zero, giving the average-ranked zero block '
        'a percentile floor of 0.065111–0.352580. Thus these percentile changes reflect '
        'a changing relative background despite unchanged zero native strength.',
        'Direction-stratified coavailability concordance was 4/8 (ρ=0.071) for '
        'astrocyte→endothelial entries and 19/24 (ρ=0.696) for the reverse direction. '
        'Native-strength concordance was 1/8 (1/4 among nonzero-both effects; ρ=−0.230) '
        'and 9/24 (9/9 among nonzero-both effects; ρ=0.817), respectively. These small, '
        'unequal subsets are descriptive and do not establish a between-direction or '
        'between-method difference. Complete scores, effects, resource attrition and '
        'tie audits are archived alongside Supplementary Figure S6.'
    ],
    'supplement_results_zh': [
        '资源审计保留原187个有向候选中的32条（17.1%）：星形胶质细胞→内皮细胞8条、'
        '内皮细胞→星形胶质细胞24条，共涉及29个不同的配体—受体定义。其余155条在'
        '原生CellChatDB中缺少完整亚基对的精确匹配；在共同全网络取交集时未再丢失'
        '目标条目。这些排除源于资源定义，并非CellChat预测失败。保留条目中31条不含'
        '原生辅助因子，1条含辅助因子，无法据此计算有意义的辅助因子分层相关性。',
        '在共同32条中，表达共可用性同向23/32（71.9%；ρ=0.565），原生强度同向'
        '10/32（31.3%；ρ=0.549），研究者计算的百分位同向16/32（50.0%；ρ=0.364），'
        '相同网络的LIANA RRA优先级同向13/32（40.6%；ρ=0.451）。原生强度在发现'
        '与外部队列分别有19与3条零效应，两队列均非零的13条中同向10/13（76.9%）。'
        'LIANA RRA分别有10与12条零效应，两队列均非零的19条中同向13/19（68.4%）。'
        '表达共可用性和百分位均无零队列效应。对这四个指标，两队列均无非零对比因'
        '10^−12容差而被置零。',
        '星形胶质细胞→内皮细胞Lamb2–Dag1、Lamb2–Itga6_Itgb1以及内皮细胞→星形胶质'
        '细胞Entpd1–Adora2b，在全部11库中的原生强度均为零，但各自的百分位对比在'
        '发现队列约为+0.014060、外部队列约为−0.024775。每库1,222个网络条目中有'
        '160–862个零强度，平均秩赋予零并列块的百分位下限为0.065111–0.352580。'
        '因此，这些百分位变化反映相对背景改变，而其原生零强度并未改变。',
        '分方向看，表达共可用性在星形胶质细胞→内皮细胞方向同向4/8（ρ=0.071），'
        '反向同向19/24（ρ=0.696）；原生强度分别同向1/8（非零效应中为1/4；'
        'ρ=−0.230）和9/24（非零效应中为9/9；ρ=0.817）。这些小且不均衡的子集'
        '仅作描述，不建立方向间或方法间差异。完整得分、效应、资源流失及并列秩审计'
        '与补充图S6一并归档。'
    ],
    'figure_s6_caption_en': (
        'Supplementary Figure S6. Controlled-resource comparison with native CellChat. '
        'Panels show cohort contrasts in (A) expression coavailability, (B) native CellChat '
        'inferred strength (net$prob), (C) analyst-derived strength percentile, and (D) '
        'LIANA unique-column RRA priority recomputed on the same 1,222-entry network. '
        'Both axes show MCAO mean minus Sham mean at the library level: GSE174574 on '
        'the horizontal axis and GSE245386 on the vertical axis. All panels retain the '
        'same 32 directed candidates. Blue circles denote the 8 astrocyte→endothelial '
        'entries and orange triangles the 24 reverse-direction entries. Hollow symbols '
        'mark the same three entries whose native strength was zero in all 11 libraries; '
        'overlapping entries may occupy identical coordinates. Annotations report '
        'Spearman ρ, concordant directions using both the complete and nonzero-effect '
        'denominators, and the number with a zero effect in either cohort. Among effects nonzero in both cohorts, '
        'native strength was concordant for 10/13 and LIANA RRA for 13/19. Hollow points '
        'in panel C illustrate percentile movement of an unchanged zero-strength block '
        'as the network background varies. The percentile is analyst-derived, not a '
        'native CellChat consensus output. These comparisons are descriptive.'),
    'figure_s6_caption_zh': (
        '补充图S6. 控制资源范围的原生CellChat比较。各面板分别展示队列对比：'
        'A，表达共可用性；B，CellChat原生推断强度（net$prob）；C，研究者计算的强度'
        '百分位；D，在相同1,222条网络背景中重算的LIANA单列RRA优先级。横轴为'
        'GSE174574、纵轴为GSE245386，均为文库层面MCAO均值减Sham均值。各面板'
        '保留相同32个有向候选。蓝色圆点表示8个星形胶质细胞→内皮细胞条目，橙色'
        '三角表示24个反向条目。空心符号在各面板均标记原生强度在全部11库中为零'
        '的相同3条；重叠条目可能位于相同坐标。图内标注Spearman ρ、以全部32条为'
        '分母及以两队列均非零效应为分母的同向计数，以及任一队列效应为零的条目数。若限定两队列效应'
        '均非零，原生强度同向10/13，LIANA RRA同向13/19。C中的空心点显示原生'
        '零强度不变时，网络背景改变仍可使零并列块的百分位位置移动。该百分位由'
        '研究者计算，并非CellChat原生共识输出。以上比较均为描述性分析。'),
    'independent_numeric_audit': {
        'path': 'work/repro_v3_cellchat_editorial_audit.json',
        'sha256': hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        'status': review['status']
    },
    'source_files': review['source_files'],
    'citation_entries': 'work/repro_v3_cellchat_references.json',
    'editorial_constraints': [
        'No general benchmark or performance-superiority claim.',
        'Native strength and analyst-derived percentile remain distinct outcomes.',
        'Resource attrition is not interpreted as failed prediction.',
        'All-zero native strengths with changing percentile effects are reported explicitly.',
        'Both complete and nonzero-effect denominators are retained.',
        'The independent implementation is not independent biological evidence.',
        'No original permutation P value is transferred to this extension.',
        'The separate full native-resource arm was not run.'
    ],
    'figure_visual_qa': {
        'status': 'pass',
        'path': 'outputs/reproducibility_v3_cellchat/figure_source_data/figureS6_independent_editorial_review.json',
        'figure_sha256': '2b79a397209af7432f6f5fa01f99dd956f49faed31790463f81ed06bafbdd109',
        'all_128_plotted_effect_rows_bitwise_equal_to_final_effects': True
    },
    'main_generators_modified': False,
}
dest = ROOT / 'repro_v3_cellchat_narrative.json'
dest.write_text(json.dumps(narrative, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({'status': narrative['status'], 'path': str(dest),
                  'main_en_words': sum(len(narrative[k].split()) for k in
                      ['main_methods_en', 'main_results_en', 'main_discussion_en'])}, indent=2))
