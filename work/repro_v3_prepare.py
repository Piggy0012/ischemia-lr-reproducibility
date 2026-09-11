"""Editorial v3: concise main manuscripts and separate complete supplements.

The v2 numerical release is immutable. Public identifiers are inserted only
from verified release metadata; unfilled fields remain conspicuous draft fields.
"""
from pathlib import Path
import json,re,shutil,hashlib

W=Path(__file__).resolve().parent
V2=W.parent/'outputs/reproducibility_v2'
OUT=W.parent/'outputs/reproducibility_v3'
OUT.mkdir(exist_ok=True)

def read(p):return p.read_text(encoding='utf-8')
def load(p):return json.loads(read(p))
def section(text, heading, stop=None):
    body=text.split(heading+'\n',1)[1]
    return body.split(stop+'\n',1)[0].strip() if stop else body.strip()
def paragraphs(v):return '\n\n'.join(v)

AUTHOR=load(W/'repro_author.json')
FRONT=load(W/'repro_manuscript_frontmatter.json')
AMBIENT=load(W/'repro_ambient_narrative.json')
THIRD=load(W/'repro_third_liana_narrative.json')
CAPS=load(W/'repro_caption_static.json')
PR='[LIANA PR #261](https://github.com/scverse/liana/pull/261)'
FIX='d4211373692e7b9c10210488ccb1efe06452b097'

ABSTRACT={
'en':('Cross-cohort expression agreement can coexist with sensitivity of ligand–receptor priority rankings to software implementation and ranking context. '
'We compared animal-level astrocyte–endothelial effects in mouse cerebral ischemia using matched candidates, dependence-preserving label permutations, controlled aggregation changes, ambient-RNA correction and a third-cohort extension. '
'Expression and native LIANA 1.10.0 priority concordance were 80.2% versus 56.1%; ranking each distinct score column once, matching a known upstream correction, raised priority-effect correlation to ρ=0.519; dynamically gated expression concordance was 79.9% versus a 50.0% joint-null median (exact P=0.005, exploratory-family Holm P=0.030), while fixed-metric and gap families did not reach adjusted significance. '
'This case study supports auditing expression effects, aggregation inputs, software versions and directional denominators separately when transferring candidate priorities across cohorts.'),
'zh':('配体—受体表达效应的跨队列一致性可与优先级排序对软件实现及排名背景的敏感性并存。'
'本研究以小鼠脑缺血的动物级星形胶质细胞—内皮细胞效应为对象，采用相同候选比较、保留依赖结构的标签置换、受控聚合改变、环境RNA校正及第三队列扩展。'
'表达与原生LIANA 1.10.0优先级同向率分别为80.2%和56.1%；每个独特评分列仅排名一次、对应已知上游修正后，优先级效应相关升至ρ=0.519；动态门槛表达同向率为79.9%，联合零分布中位数为50.0%（精确P=0.005，探索家族Holm P=0.030），固定指标与差值家族均未达校正显著。'
'本案例支持在跨队列转移候选优先级时，分别审计表达效应、聚合输入、软件版本及方向性分母。')}

METHODS={
'en':[
('Data and cell selection',[
'GSE174574 contributed three sham and three acute middle cerebral artery occlusion (MCAO) animal libraries; GSE245386 contributed three sham and two newly generated wild-type MCAO libraries, excluding knockout samples and reused data. Samples were reconciled with GEO and the source papers, including the published correction to GSE245386. The animal represented by each original library was the comparison unit. [@Zheng2022][@Ruan2023][@RuanCorrection2025]',
'Matrices were processed separately by library. Quality control retained 200–6,000 detected genes, at least 500 UMIs and at most 20% mitochondrial UMIs. A fixed marker rule assigned broad lineages and retained unassigned cells. In the original cohorts, an alternative selection intersected the rule labels with an agreeing CellTypist whole-mouse-brain label, reference score ≥0.5 and no Scrublet-predicted doublet. Reference scores were classifier outputs rather than calibrated identity probabilities. Full marker rules, models and diagnostic parameters are in the Reproducibility appendix. [@Dominguez2022][@Yao2021][@Yao2023][@Wolock2019]',
'GSE332910 supplied six ipsilateral striatal single-nucleus RNA libraries after sham surgery or acute MCAO. The source supplementary table reports three biological replicates and three libraries per group, without specifying mouse or pool composition within each library. The unchanged marker rule was applied, and this extension was analyzed descriptively outside the original permutation families. [@Zhang2026]']),
('Expression scores and relative priorities',[
'Using a fixed mouseconsensus resource, we evaluated astrocyte-to-endothelial candidates and the reverse direction. Per-animal expression coavailability was the square root of the product of ligand and receptor mean log1p(CP10k); complexes used the minimum constituent value. Disease effects were equal-animal MCAO means minus sham means. The dynamic expression gate required all subunits to be detected in at least 10% of the relevant cells in at least two animals in either condition, assessed separately for ligand and receptor.',
'LIANA 1.10.0 was run per library on all eligible assigned lineages before the two target directions were extracted. CellPhoneDB, Connectome, log2FC, NATMI and SingleCellSignalR were method interfaces within this workflow. Coavailability, lr_means, expr_prod and lrscore were grouped as related expression summaries; Connectome and NATMI shared expr_prod. Relative priority was 1−magnitude_rank, where magnitude_rank is a bounded robust rank aggregation score. The interfaces therefore supplied related score representations rather than independent replications. [@Dimitrov2022][@Dimitrov2024][@Kolde2012]',
'The main comparison fixed candidates with finite values for every metric in all original animal libraries, without a treatment-dependent additional gate. Directional concordance used both all candidates and, separately, those with nonzero effects in both cohorts as denominators; zero–zero pairs were not concordant. Spearman correlation included all candidates. Missing scores were not replaced by zero. Third-cohort comparisons used the same complete candidates across all original and third-cohort libraries, with rankings calculated before target restriction.']),
('Label randomization and multiplicity',[
'Whole-animal score vectors were held intact while condition labels were exhaustively reassigned within each original study. The 20 discovery and 10 validation allocations yielded 200 joint configurations, including the observed labels, preserving shared-gene dependence. Each configuration recomputed concordance and correlation; the dynamic-gate analysis also reapplied expression eligibility. One-sided exact P values were upper-tail fractions of these configurations. Null percentiles describe a random-label reference interval rather than uncertainty intervals for observed concordance.',
'Holm correction was applied separately to family A (30 fixed-metric tests), family B (24 expression-minus-native-priority gap tests) and family C (six exploratory dynamic-gate tests). The families were specified in a revision-stage plan after earlier results had been inspected and were not preregistered. The statistics share inputs and are dependent; Holm controls family-wise error under arbitrary dependence and was retained for the original six-test family even where statistics were redundant. Equal P values alone do not establish dependence. The gap null concerns treatment-label association, not algorithm equivalence.',
'A supplementary conditional diagnostic held one cohort at its observed labels while reallocating the other. Its resolution was limited by the separate allocation spaces, including a minimum validation P of 0.10. The finer joint null did not increase the number of biological replicates or establish disease association separately in both cohorts.']),
('Controlled aggregation and ranking context',[
'A deterministic example and complete per-library networks were used to compare the installed aggregation with a diagnostic that ranked each distinct magnitude-score column once. We independently executed the aggregation function from the fixed merge commit of '+PR+' on the same saved full-network scores. This isolated the known aggregation correction; it was not a full newer-release pipeline rerun. Resource membership, detection threshold and upstream scores were held unchanged.',
'We then restricted aggregation to the full-network edge intersection shared by all original libraries before extracting the unchanged target set. A separate fixed-target percentile transform evaluated a simpler ranking representation. Network membership, saturation, zero effects and both directional denominators were reported. These additions were descriptive; permutation P values from the original pipeline were not transferred to them.']),
('Ambient RNA sensitivity',[
'DecontX was fitted separately to all QC-passing cells in each original library using existing whole-brain reference labels as clusters and no measured empty-droplet background. Numerical convergence was checked, and the single unconverged default fit was extended under unchanged inputs and settings. All final selected fits reached the specified threshold. Model settings and the complete default-versus-extended comparison are in the Reproducibility appendix. [@Yang2020]',
'Corrected fractional counts were normalized by corrected full-gene totals while cell selections remained fixed. Detection greater than zero was the main corrected rule; estimated count ≥1 was a separate eligibility sensitivity. Original and corrected scores were also compared on identical complete candidates. These estimates were not rounded into observed UMIs or supplied to DESeq2. Count subtraction and retained eligibility remain conditional on model and cluster assumptions.'])],
'zh':[
('数据与细胞选择',[
'GSE174574纳入3个假手术及3个急性大脑中动脉闭塞（MCAO）动物文库；GSE245386纳入3个假手术及2个新生成的野生型MCAO文库，排除敲除样本与重复使用的数据。依据GEO、原论文及GSE245386已发表更正核对样本，以每个原始文库代表的动物为比较单位。[@Zheng2022][@Ruan2023][@RuanCorrection2025]',
'矩阵逐文库处理，保留检出基因200–6,000、UMI至少500且线粒体UMI不超过20%的细胞。固定标记规则分配大谱系并保留未分配细胞。原两个队列另设参考支持选择，要求规则标签与CellTypist全鼠脑标签一致、参考得分至少0.5且未被Scrublet标记为双细胞。模型得分不解释为校准后的身份概率；完整标记规则、参考模型与参数见复现附录。[@Dominguez2022][@Yao2021][@Yao2023][@Wolock2019]',
'GSE332910提供假手术或急性MCAO后的6个同侧纹状体单核RNA文库；来源补充表报告每组3个生物学重复及3个文库，每文库单鼠或混样构成未说明。沿用相同标记规则，以描述性扩展分析，不纳入原两队列的置换家族。[@Zhang2026]']),
('表达评分与相对优先级',[
'使用固定mouseconsensus资源，评价星形胶质细胞至内皮细胞及其反向候选。逐动物表达共可用性为配体与受体平均log1p(CP10k)之积的平方根，复合物取最小亚基值；疾病效应为MCAO动物均值减Sham动物均值。动态表达资格要求所有亚基在至少一组的至少2个动物中、相应细胞检出比例达到10%，配体与受体分别判断。',
'LIANA 1.10.0逐文库在全部合格大谱系上计算后提取目标两方向。CellPhoneDB、Connectome、log2FC、NATMI及SingleCellSignalR是同一流程中的方法接口。共可用性、lr_means、expr_prod与lrscore归为相关的表达汇总，其中Connectome与NATMI共享expr_prod。相对优先级为1减magnitude_rank，后者是有界的稳健秩聚合分数。因此各接口提供相关的评分表示，不构成独立重复验证。[@Dimitrov2022][@Dimitrov2024][@Kolde2012]',
'主比较固定全部原始动物、全部指标均完整的共同候选，不额外使用依赖处理标签的门槛。同向率分别以全部候选及两队列效应均非零的候选为分母，零—零不记为同向；Spearman相关包含全部候选，缺失评分不补零。第三队列使用全部原始及新增文库共同完整的相同目标候选，优先级均在限制目标前完成排名。']),
('标签随机化与多重检验',[
'固定每个动物完整的评分向量，在原两个研究内穷举处理标签，发现队列20种、验证队列10种配置形成含观测标签在内的200种联合配置，以保留共享基因依赖。每种配置重算同向率与相关；动态门槛分析同时重算表达资格。单侧精确P为配置统计量的上尾比例，零分布分位范围是随机标签参照区间。',
'分别在A家族30项固定指标、B家族24项表达减原生优先级差值及C家族6项探索性动态门槛统计内使用Holm校正。家族在查看旧结果后的修订计划中规定，并非前瞻注册。统计量共享输入而具有依赖，Holm可在任意依赖下控制家族错误率；即使部分统计冗余，仍保留原6项家族。相同P值本身不能证明相关。差值零假设针对处理标签关联，不等于算法等效性检验。',
'补充条件诊断固定一队列观测标签，仅重分配另一队列；其精度受各自配置数限制，验证队列最低P仍为0.10。联合配置更细不增加生物学重复，也不独立确立两个队列各自的疾病关联。']),
('受控聚合与排名背景',[
'通过确定性示例和逐文库完整网络，对比已安装的聚合实现与每个独特幅度评分列只排名一次的诊断。另在相同完整网络评分上独立执行'+PR+'固定合并提交的聚合函数，隔离已知聚合修正；这一检查不作为完整新版流程重跑。资源成员、检测门槛和上游评分保持不变。',
'随后在全部原始文库共同的全网络边交集内聚合，再提取相同目标集合；另以固定目标百分位转换评价更简单的排名表示。报告网络成员、聚合饱和、零效应及两种方向分母。这些追加比较均为描述性，原生流程的置换P不移用于追加结果。']),
('环境RNA敏感性',[
'DecontX逐原始文库拟合全部质控合格细胞，以既有全脑参考标签分群，不提供实测空液滴背景。逐库检查数值收敛，唯一默认未收敛文库在相同输入与设置下延长拟合；最终选用的全部结果达到预设阈值。参数与默认及延长拟合的完整比较见复现附录。[@Yang2020]',
'保持细胞选择不变，分数校正计数按全部基因校正总量归一化。大于零为主检测规则，估计计数至少1作为独立的资格敏感性；并在相同完整候选上比较原始及校正评分。估计值不取整为观测UMI，也不输入DESeq2；扣减量和资格保留均受模型及分群假设约束。'])]}

CHECKLIST={
'en':[
'Preserve the biological unit: retain per-animal expression scores, condition labels and each library’s mouse or pool composition.',
'Record cell selection: publish marker panels, reference-label mappings, thresholds, doublet decisions and retained cell counts.',
'Export candidate eligibility: retain per-subunit detection, missingness and the candidate set used in each comparison and permutation.',
'Pin software and resources: record package versions, source commits, interaction resources and explicit parameters; rerun version-sensitive checks.',
'Preserve aggregation inputs: retain the complete network, distinct magnitude columns and their ranking directions before consensus calculation.',
'Report ranking context: list network membership, tie rules and saturation counts, and distinguish changing backgrounds from fixed intersections.',
'State directional denominators: report all-candidate and nonzero-effect agreement, correlation and zero-effect counts alongside all prespecified test families.'
],
'zh':[
'保留生物学单位：输出逐动物表达分数、处理标签及每文库单鼠或混样构成。',
'记录细胞选择：提供标记面板、参考标签映射、阈值、双细胞判定及最终细胞数。',
'导出候选资格：保留逐亚基检出比例、缺失状态，以及每个比较和每次置换使用的候选集合。',
'固定软件与资源：记录包版本、源码提交、相互作用资源及显式参数，重跑受版本影响的检查。',
'保留聚合输入：输出共识计算前的完整网络、独特幅度评分列及其排名方向。',
'报告排名背景：提供网络成员、并列处理与饱和数量，区分变化背景和固定交集。',
'列明方向分母：同时报告全部候选及非零效应同向率、相关、零效应数，并完整列出预先规定的检验家族。']}

def statements(lang):
    d=load(W/'repro_v3_author_declarations.json')
    rec=load(OUT/'archive_record.json') if (OUT/'archive_record.json').exists() else {}
    public=rec.get('public_verified') is True
    if public:
        assert rec['repository_url'].startswith('https://github.com/Piggy0012/')
        assert re.fullmatch(r'10\.5281/zenodo\.\d+',rec['doi'])
        assert re.fullmatch(r'[a-f0-9]{40}',rec['analysis_commit'])
        release=f'[{rec["repository_url"]}]({rec["repository_url"]}); [version DOI {rec["doi"]}](https://doi.org/{rec["doi"]}); analysis commit `{rec["analysis_commit"]}`.'
    elif rec.get('github_verified') is True:
        assert rec['repository_url'].startswith('https://github.com/Piggy0012/')
        assert re.fullmatch(r'[a-f0-9]{40}',rec['analysis_commit'])
        release=f'[{rec["repository_url"]}]({rec["repository_url"]}); [VERSION_SPECIFIC_ZENODO_DOI]; analysis commit `{rec["analysis_commit"]}`.'
    else:
        release='[PUBLIC_REPOSITORY_URL]; [VERSION_SPECIFIC_ZENODO_DOI]; [ANALYSIS_COMMIT].'
    if lang=='en':
        code='## Code availability\n\n'+release+'\n\nThe release comprises analysis scripts, versioned environment files, cell-summary counts, complete result tables and figure source data. The Reproducibility appendix and repository run instructions distinguish the full LIANA 1.10.0 workflow from replay of the aggregation function at upstream commit `'+FIX+'` on saved complete networks; the latter does not represent a full fixed-commit pipeline rerun. The file index maps each named input and output to its repository path.\n\n'
        funding=d.get('funding_statement_en') or '[FUNDING_AGENCY_AND_GRANT_NUMBER_TO_BE_CONFIRMED]'
        return code+'## Funding\n\n'+funding+'\n\n## Competing interests\n\nThe author declares no competing interests.\n\n## Author contributions\n\nSihuan Zhu: conceptualization, research direction, critical revision, approval of the final manuscript and responsibility for the work. AI-assisted computational and writing activities are disclosed below.\n\n## Declaration of AI use\n\nOpenAI Codex was used to assist study planning, public-data retrieval, code development and execution, statistical and figure checks, and drafting, translation and revision of the manuscript. Numerical results were calculated by the archived analysis code from the stated source data. Sihuan Zhu reviewed and approved the manuscript and takes responsibility for its content.\n\n## Acknowledgments\n\nThe author acknowledges the investigators who released the source datasets and the developers of the analytical software, including the LIANA maintainers for the aggregation correction documented in '+PR+'.\n\n'
    funding=d.get('funding_statement_zh') or '[资助方与项目编号待确认]'
    return '## 代码获取\n\n'+release+'\n\n发布内容包括分析脚本、版本化环境文件、细胞汇总计数、全部结果表及图源数据。复现附录与仓库运行说明区分LIANA 1.10.0完整流程，以及在已保存完整网络上重放上游提交`'+FIX+'`的聚合函数；后者不作为完整修复提交流程重跑。文件索引将文中引用的输入和输出对应到仓库路径。\n\n## 经费支持\n\n'+funding+'\n\n## 利益冲突\n\n作者声明不存在利益冲突。\n\n## 作者贡献\n\n朱四欢：研究构思、研究方向确定、重要内容修订、最终稿批准及研究责任承担。AI辅助的计算与写作活动如下披露。\n\n## AI使用声明\n\n本研究使用OpenAI Codex辅助研究规划、公开数据获取、代码开发与执行、统计及图表检查，以及稿件起草、翻译和修订。数值结果由归档分析代码依据所列来源数据计算。朱四欢审阅并批准稿件，对其内容承担责任。\n\n## 致谢\n\n感谢公开来源数据的研究者及分析软件开发者，并感谢LIANA维护者在'+PR+'中公开聚合修正。\n\n'

def main():
    ccpath=W/'repro_v3_cellchat_narrative.json'
    cc=load(ccpath) if ccpath.exists() else None
    if cc is not None: assert cc['status']=='complete'
    for folder in ['figures','figure_source_data']:
        shutil.copytree(V2/folder,OUT/folder,dirs_exist_ok=True)
    if cc:
        for item in cc['source_files']:
            assert hashlib.sha256((W.parent/item['path']).read_bytes()).hexdigest()==item['sha256'],item['path']
        figure_root=W.parent/'outputs/reproducibility_v3_cellchat/figures'
        basename='FigureS6_cellchat_controlled_extension'
        assert cc['figure_visual_qa']['status']=='pass'
        assert hashlib.sha256((figure_root/(basename+'.png')).read_bytes()).hexdigest()==cc['figure_visual_qa']['figure_sha256']
        for ext in ['png','pdf','svg']:
            shutil.copy2(figure_root/(basename+'.'+ext),OUT/'figures'/(basename+'.'+ext))
    for lang in ['en','zh']:
        src=read(W/f'repro_manuscript_{lang}_source.md')
        part=load(W/('repro_manuscript_en_sections.json' if lang=='en' else 'repro_main_sections.json'))
        en=lang=='en'
        headers=['Abstract','Introduction','Methods','Results','Discussion','Conclusion','Data availability','Figure legends','References'] if en else ['摘要','引言','方法','结果','讨论','结论','数据获取','图注','参考文献']
        intro=part['introduction'][:2]+[part['introduction'][-1]]
        result=section(src,'## '+headers[3],'## '+headers[4])
        # Tables are already source-validated by the v2 assembler.
        table_anchor='Table 1.' if en else '表1'
        result,tables_text=result.split(table_anchor,1)
        tables_text=table_anchor+tables_text
        # Auxiliary datasets appear once in the body, exclusively as supplements.
        if en:
            bio=result.split('Sorted-data summaries',1)[1].strip()
            bio='Sorted-data summaries '+bio
            result=result.split('Sorted-data summaries',1)[0].strip()
            result=result.replace('### A third cohort and sorted expression define the biological scope','### Transfer to a third cohort')
            result=result.replace('All six dynamically gated statistics had exact P=0.005 and Holm-adjusted P=0.030 within exploratory family C (Figure 2D).','All six dynamically gated statistics had exact P=0.005 and Holm-adjusted P=0.030 within exploratory family C (Figure 2D). For each statistic, the observed allocation was the unique upper-tail maximum among all 200 allocations, so P=1/200 was the exact resolution limit. Within each selection, the two concordance fractions were identical across allocations because all effects were nonzero; the six-test correction was retained.')
            astart=result.index('DecontX fitting was completed')
            aend=result.index('In primary-rule target cells',astart)
            result=result[:astart]+('All eleven final DecontX fits met the specified numerical threshold. The pooled cell-level median estimated contamination was 6.13%; these cell-level summaries were descriptive. Complete fit diagnostics and the single extended fit are reported in the Reproducibility appendix.\n\n')+result[aend:]
            result+='\n\nBarrier-related transcription, paired-pool sorted expression and spatial descriptions are provided solely as supplementary biological context (Supplementary analyses and Figures S3–S5).'
        else:
            bio='### 分选效应'+result.split('### 分选效应',1)[1]
            result=result.split('### 分选效应',1)[0].strip()
            result=result.replace('预先划定的6项动态门槛探索统计经Holm校正均为0.030（图2D）。','修订计划规定的6项动态门槛探索统计经Holm校正均为0.030（图2D）。每项统计的观测配置均为全部200种配置中唯一的上尾最大值，故P=1/200达到精确分辨率下限。同一细胞选择中的两种同向率在全部配置下完全相同，因为所有效应均非零；仍保留原6项家族校正。')
            astart=result.index('完成了 11 个文库')
            aend=result.index('原规则指定的目标细胞中',astart)
            result=result[:astart]+'全部11份最终DecontX拟合达到预设数值阈值，合并细胞的估计污染中位数为6.13%，该细胞层级汇总仅作描述。完整拟合诊断及唯一延长拟合的记录见复现附录。\n\n'+result[aend:]
            result+='\n\n屏障相关转录、配对样本池分选表达及空间描述仅作为生物学背景，全部置于补充分析与图S3–S5。'
        # Keep the first three cohort rows in the main table.
        tables_text='\n'.join(line for line in tables_text.splitlines() if not line.startswith('| GSE163752') and not line.startswith('| GSE233814'))
        if en:
            disc=section(src,'## Discussion','## Conclusion').split('\n\n')
            disc=[p for p in disc if not p.startswith('The practical implication')]
            disc=[p for p in disc if not p.startswith('Ambient-RNA sensitivity requires')]
            disc=[p.replace('while sorted ipsilateral–contralateral contrasts and spatial descriptions interrogate different biological comparisons. ','') for p in disc]
            disc=[p.replace('Cell identity and contamination sensitivity require their own audits, Functional ligand release', 'Cell identity and contamination sensitivity require their own audits. Functional ligand release') for p in disc]
            disc.insert(-1,'Ambient-count estimation adds a distinct source of sensitivity. The corrected analysis changed correlations on matched candidates, while an estimated-count threshold of at least one selected a much smaller subset. This shows why count estimation and candidate eligibility should be retained as separate outputs: an apparent gain in agreement after thresholding can reflect a change in the assessed candidates.')
            checklist_title='### A practical reporting checklist'
        else:
            disc=section(src,'## 讨论','## 结论').split('\n\n')
            disc=[p for p in disc if not p.startswith('环境RNA敏感性分析需要')]
            disc=[p.replace('分选样本采用池化和对侧对照，空间数据没有独立区域或动物层面的验证，且','') for p in disc]
            disc.insert(-1,'环境计数估计带来另一层敏感性：相同候选上的相关随校正而变，而估计计数至少1的门槛保留了明显更小的集合。因此，应分别记录计数估计与候选资格，才能区分评分本身的改变和筛选对象改变所造成的同向率变化。')
            checklist_title='### 可直接使用的报告清单'
        if cc:
            result+='\n\n### '+('A controlled independent-framework extension' if en else '受控的独立框架扩展')+'\n\n'+cc['main_results_'+lang]
            disc.insert(-1,cc['main_discussion_'+lang])
        discussion=paragraphs(disc)+'\n\n'+checklist_title+'\n\n'+'\n\n'.join(f'{i}. {s}' for i,s in enumerate(CHECKLIST[lang],1))
        head=src.split('## '+headers[0],1)[0]
        abstract=ABSTRACT[lang]
        if cc:
            abstract=abstract.replace('ambient-RNA correction and a third-cohort extension','ambient-RNA correction, a third-cohort extension and a controlled native CellChat comparison').replace('及第三队列扩展','、第三队列扩展及受控的原生CellChat比较')
        maintext=head+'## '+headers[0]+'\n\n'+abstract+'\n\n'+('Keywords: ' if en else '关键词：')+FRONT['keywords_'+lang]+'\n\n'
        maintext+='## '+headers[1]+'\n\n'+paragraphs(intro)+'\n\n## '+headers[2]+'\n\n'
        for h,body in METHODS[lang]:maintext+='### '+h+'\n\n'+paragraphs(body)+'\n\n'
        if cc:
            cc_methods=cc['main_methods_'+lang].replace('新增队列对比仅作描述性分析','上述新增比较仅作描述性分析')
            maintext+='### '+('Independent CellChat comparison' if en else '独立CellChat比较')+'\n\n'+cc_methods+'\n\n'
        maintext+='## '+headers[3]+'\n\n'+result+'\n\n'+tables_text+'\n\n## '+headers[4]+'\n\n'+discussion+'\n\n'
        maintext=maintext.replace('Equal P values alone do not establish dependence. ', '').replace('相同P值本身不能证明相关。', '')
        maintext+='## '+headers[5]+'\n\n'+FRONT['conclusion_'+lang]+'\n\n'
        links=', '.join(f'[{a}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a})' for a in ['GSE174574','GSE245386','GSE332910','GSE163752','GSE233814'])
        maintext+='## '+headers[6]+'\n\n'+('Public source data are available from GEO: ' if en else '原始公开数据可由GEO获取：')+links+'.\n\n'+statements(lang)
        maintext+=('## Ethics statement\n\nThis study reanalyzed public animal data and performed no new animal or human experiments. Source approvals are reported in the original studies.' if en else '## 伦理说明\n\n本研究仅重分析公开动物数据，未开展新的动物或人体实验。原始实验审批见来源论文。')+'\n\n'
        maintext+='## '+headers[7]+'\n\n'
        legends={**CAPS[lang],'Figure4':AMBIENT['caption_'+lang],'Figure5':THIRD['caption_'+lang]}
        for key in ['Figure1','Figure2','Figure3','Figure4','Figure5']:maintext+=legends[key]+'\n\n'
        maintext+='## '+headers[8]+'\n\n<!-- REFERENCES_GENERATED -->\n'
        appendix='# '+('Reproducibility appendix and supplementary analyses' if en else '复现附录与补充分析')+'\n\nSihuan Zhu\n\n'
        appendix+=('This supplement contains complete operational methods, numerical provenance and auxiliary biological analyses for the companion manuscript. All reported null families and diagnostic outputs remain separately identified.' if en else '本补充材料提供主稿的完整操作方法、数值溯源及辅助生物学分析。各置换家族及诊断输出保持独立标识。')+'\n\n'
        appendix+='## '+('Reproducibility appendix' if en else '复现附录')+'\n\n'
        appendix+=section(src,'## '+headers[2],'## '+headers[3])+'\n\n'
        if en:
            appendix=appendix.replace('GSE233814 contributed five section records to a supplementary expression figure using verified spatial coordinates.', 'GSE233814 contributed five section records to archived descriptive analyses; the supplementary expression figure displays the Control and Day 1 sections using verified spatial coordinates.')
        appendix+='### '+('Complete fit and numerical diagnostics' if en else '完整拟合与数值诊断')+'\n\n'+AMBIENT['results_'+lang].split('\n\n')[1]+'\n\n'
        if cc:
            appendix+='### '+('Independent CellChat methods and provenance' if en else '独立CellChat方法与溯源')+'\n\n'+paragraphs(cc['supplement_methods_'+lang])+'\n\n'
        appendix+='### '+('Sources and release layout' if en else '来源与发布结构')+'\n\n'
        appendix+=('Repository paths and checksums are listed in FILE_INDEX.tsv; figure inputs are mapped in FIGURE_SOURCE_INDEX.tsv. ENVIRONMENT.md distinguishes the tested analysis runtime from restoration recipes. RUNNING.md describes cached reconstruction, raw inference and fixed-commit aggregation replay. Public version identifiers are provided in the companion Code availability section.' if en else '仓库路径与哈希列于FILE_INDEX.tsv，图源数据由FIGURE_SOURCE_INDEX.tsv映射。ENVIRONMENT.md区分实际验证环境和恢复方案；RUNNING.md区分缓存重建、原始推断与固定提交聚合函数重放。公开版本标识列于主稿代码获取段。')+'\n\n'
        if cc:
            appendix+='## '+('Controlled independent-framework results' if en else '受控独立框架结果')+'\n\n'+paragraphs(cc['supplement_results_'+lang])+'\n\n'
        appendix+='## '+('Supplementary biological analyses' if en else '补充生物学分析')+'\n\n'+bio+'\n\n'
        appendix+='## '+('Supplementary figure legends' if en else '补充图注')+'\n\n'
        for key in ['FigureS1','FigureS2','FigureS3','FigureS4','FigureS5']:appendix+=legends[key]+'\n\n'
        if cc:appendix+=cc['figure_s6_caption_'+lang]+'\n\n'
        appendix+='## '+headers[8]+'\n\n<!-- REFERENCES_GENERATED -->\n'
        for name,txt in [('manuscript',maintext),('supplement',appendix)]:
            (W/f'repro_v3_{name}_{lang}_source.md').write_text(txt,encoding='utf-8')
        main_methods=section(maintext,'## '+headers[2],'## '+headers[3])
        for token in ['SHA256','float32','CSR','iterLogLik','20260911','marker_panel.json']:
            assert token not in main_methods,(lang,token)
        assert 'Code availability' in maintext if en else '代码获取' in maintext
    rec=load(OUT/'archive_record.json') if (OUT/'archive_record.json').exists() else {}
    public=rec.get('public_verified') is True
    fields=[] if public else (['VERSION_SPECIFIC_ZENODO_DOI'] if rec.get('github_verified') is True else ['PUBLIC_REPOSITORY_URL','VERSION_SPECIFIC_ZENODO_DOI','ANALYSIS_COMMIT'])
    (OUT/'manuscript_readiness.json').write_text(json.dumps({'submission_ready':False,'public_archive_verified':public,'author_approval_supplied':True,'competing_interests_supplied':True,'funding_details_supplied':bool(load(W/'repro_v3_author_declarations.json').get('funding_statement_en')),'draft_fields':fields,'note':'Final rendered-document and identifier verification are recorded separately; this source assembly does not establish submission readiness.'},indent=2),encoding='utf-8')
    print('Prepared v3 main and supplement sources without modifying v2 scientific results')

if __name__=='__main__':main()
