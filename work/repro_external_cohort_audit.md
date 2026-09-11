# 外部队列扩展与方法学重叠审计

核验日期：2026-09-10。原始稿件的两个单细胞队列及既有排除项沿用 `work/dataset_audit.json`；本次重新检索 GEO，不把新发表的二次分析视为新样本。检索与机器可读清单见 `repro_external_cohort_audit.json`。

## 可实施选择

**GSE332910 是本次找到的最优急性第三队列候选。**其于2026-09-09公开，原论文为 Zhang 等的 Nrsn1–Smarcc1 研究。GEO明确列有3个24 h RNA文库和3个Sham RNA文库；原文确认60 min MCAO再灌注、同侧纹状体、C57BL/6小鼠，并识别星形胶质及内皮细胞。其余0 h、7 d和14 d，以及对应ATAC记录，不计入本次急性RNA比较。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332910)，[原论文](https://doi.org/10.1002/advs.77547)。

**补充材料已明确作者报告的生物学重复，个体鼠映射仍待补全。**2026-09-11经PMC新的公开Cloud通道取得XLSX与DOCX附件，MD5与官方metadata逐一核对。Table S1的Sham和24 h两列均明确写“Biological replicates=3”“Number of libraries=3”。因此准确措辞是“作者报告每组3个生物学重复、3个文库”；不能据此单独确认逐文库各为一只鼠还是独立混样。正文、补充说明及对应6条BioSample均未提供动物ID、pool大小或文库级性别。该队列可作为新的急性生物学重复扩展，同时保留单核纹状体与原scRNA队列的采样差别及GSE245386原来的n=2事实。

最初网页附件返回HTML检查，旧oa.fcgi接口返回404。进一步核验发现PMC在2026-08-26完成旧OA API/FTP分发下线，新公开S3 metadata提供了可匿名取回附件的正式地址，并非绕过访问限制。[官方Cloud说明](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/)，[Table S1原始XLSX](https://pmc-oa-opendata.s3.amazonaws.com/PMC13542395.1/ADVS-9999-e77547-s001.xlsx)。来源、校验和与完整提取表格见 `work/repro_cohort_metadata/GSE332910_supplement_download_audit.json` 和 `ADVS-9999-e77547-s001.json`。

| 组别 | 文库 | GEO提交压缩文件大小 |
| --- | --- | ---: |
| 24 h | GSM9755192 | 137,766,625 bytes |
| 24 h | GSM9755194 | 101,499,743 bytes |
| 24 h | GSM9755196 | 93,333,915 bytes |
| Sham | GSM9755210 | 88,263,041 bytes |
| Sham | GSM9755212 | 101,326,327 bytes |
| Sham | GSM9755214 | 39,760,025 bytes |

六个文件已完整下载，合计561,949,676 bytes，文件大小及SHA-256均已记录。完整URL见 `repro_cohort_metadata/GSE332910_download_manifest.json`。压缩包内部为 `barcodes.tsv.gz`、`features.tsv.gz` 和 `matrix.mtx.gz`；PISA产生整数UMI稀疏矩阵，features只有一列提交基因符号，不能误读为Ensembl两列格式。barcode为CELL编号，必须以GSM区分文库。矩阵按基因行排序，不能套用旧10x数据的“按细胞列递增”断言。

第一个24 h文库头为30,125基因×25,095细胞、47,662,066个非零条目；第一个Sham为29,877基因×18,445细胞、30,971,985个非零条目。它们是提交矩阵尺寸，不是质控后细胞数。标准化文件仅将提交符号复制成identifier与symbol两列，保留来源含义，不伪造Ensembl映射。初步提取见 `work/repro_third_cohort.py`；该脚本复用原始 `stream_annotate.py` 的QC与谱系规则，并另存下载清单、校验和、映射说明及每文库细胞汇总。

区域、核提取和BGI平台均与原来两组全半球scRNA不同。因此，加入该队列主要评价跨研究及采样技术的可迁移程度；若结果不同，不能归结为生物学时窗不同，也不能直接证明某种工具错误。

## 已完成的第三队列分析

2026-09-11完成六文库的下载、相同参数QC及谱系打分。119,278个提交核中116,308个通过QC，原规则指定6,430个星形胶质核和1,112个内皮核，每文库至少75个内皮核；77,582个通过QC的核（66.70%）仍未获谱系标签。高未分配率提示原marker规则向核转录组转移有限，当前标签不能等同于作者注释或独立参考鉴定。逐文库表见 `work/repro_third_cohort/library_qc_summary.tsv`，完整细胞分数、阈值及来源见该目录README和JSON。

脚本 `work/repro_third_coavailability.py` 分块处理全部6个gene-major矩阵，验证基因-细胞坐标唯一性、非零数及已分配细胞总UMI守恒，生成全部15类每文库原始pseudobulk、平均log1p(CP10k)和检出比例。以6个文库均有的28,483个提交符号建立覆盖集合，不将未提交的符号补零。双方向共有7,582条资源覆盖候选，5%、10%和20%资格阈值分别保留1,417、828和368条。

按原10%资格规则，在与GSE174574共同合格的223条中141条同向（63.23%，Spearman ρ=0.560）；与GSE245386共同合格的445条中267条同向（60.00%，ρ=0.547）。这是每次使用不同共同候选集合的描述性比较，不能将差异特定归因为平台或区域；独立动物数未核实，不计算动物层P，也不将第三队列加入原两队列200标签置换。完整分数、全部候选组差、资格和交集行表在 `outputs/reproducibility_v2/tables/third_*`。

五个固定示例中，Timp3–Kdr、Ptn–Ptprz1及Plat–Lrp1通过10%和20%资格，组差分别为+0.04331、−0.00406和+0.27817；Ptn–Ptprz1接近零，不宜仅凭符号赋予生物学方向。Spp1–Itga5_Itgb1与Col4a1–Itga3_Itgb1虽呈正组差（+0.07911和+0.04393），均未达到10%资格，仅在5%阈值通过。Spp1在星形核的缺血组检出率为5.41%、11.70%和2.40%，只有1个库达到10%。这提供表达资格的外部边界，不是Spp1通讯机制验证。

## 完整网络LIANA及共同目标集合比较

第三队列6个LIANA完整网络均已运行完成，合计224,004条网络行及4,194条合格目标行。原生rank逐库从完整网络准确重建，并在同一网络上计算每个不同分数列仅排名一次后实施RRA的unique-column诊断；原始稀疏矩阵归一化后所有存储值均为正。与原两队列采用相同99条、全部17文库及全部相关指标均完整的目标候选比较，原发现与原外部之间的custom/native/diagnostic同向数为80/65/70，ρ为0.767/0.350/0.571。第三队列与原发现之间依次为53/48/52，ρ为0.551/0.219/0.460；与原外部之间为58/48/48，ρ为0.576/0.287/0.392，三对比较分母均为99。排序始终在逐文库完整背景网络上计算，99为比较目标子集。

这支持将方向一致率与效应相关性分开：诊断提高了第三队列优先级效应ρ，同向比例却仅小幅变化或不变；表达方向的一致性本身也在跨脑区/单核技术转移时下降。固定目标集合控制了比较候选覆盖的差别，但未分离平台、组织区域、处理或背景网络的具体作用。全部结果均为描述性，第三队列仅有raw/primary版本。完整数据在 `outputs/reproducibility_v2/tables/third_rank_diagnostic_*`，图5与双语文字见相应Figure5文件及 `work/repro_third_liana_narrative.json`。

## 筛选范围与主要排除证据

GEO GDS检索式为 `(stroke OR cerebral ischemia OR MCAO) AND (single cell OR single nucleus) AND Mus musculus[Organism] AND gse[Entry Type]`，返回145条（包含bulk、空间及非脑命中）。保存全部检索结果，并对21个新候选取回完整SOFT，结合既有审计优先检查能够增加独立对照/缺血文库且同时采到两目标细胞的候选。这是有明确范围的更新检索，不宣称PRISMA系统综述或穷尽全部可能数据。

| 数据 | 可核实结构 | 不作为更强24 h替代的原因 |
| --- | --- | --- |
| GSE225948 | young D2 brain有4个Sham和4个stroke重复文库 | 取样为CD45免疫细胞及CD45−/Ly6c高内皮；未系统采样astro。GEO提取方法明确每次实验混合4–5个半球。其余blood文库不能转成brain重复。 |
| GSE267240 | D3，雌雄各一Sham/一stroke文库；每文库4只鼠混样 | 特异分选microglia，缺两目标细胞且每性别条件无文库重复。 |
| GSE319237 / GSE313837 / GSE268505 | sorted microglia或细胞系/组织bulk RNA | 不是本研究的astro–EC单细胞比较。 |
| GSE279665 / GSE279666 | 前者6条brain_1–6未在GEO映射条件；后者Healthy、Sham、D3、D14各一snRNA文库 | 原文研究D3/D14；无24 h，不能从总数6推定3对3。 |
| GSE289791 | WT cortex及WT ipsi/contra各一文库，另有KO对应项 | 无独立重复Sham，左右半球不是额外动物。 |
| GSE247102 | D3白质WT ipsi 2、contra 2，另有KO | 无Sham且病例文库数仍为2。 |
| GSE303092 / GSE290194 | D14 ILC2或PLX5622干预，部分whole-brain 2文库 | 缺同批次Sham，不能把干预对照改称未缺血。 |
| GSE300442 | 分选内皮；不同半球/处理各一文库 | 缺astro及重复。 |
| GSE295882 | ATAD3A flox 3文库与astro KO 4文库 | 基因型比较，未提交可追踪的WT Sham vs MCAO两组。 |
| GSE300564 / GSE310324 / GSE319238 / GSE335511 | 每治疗/区域/饮食条件一文库，部分为免疫富集 | 没有符合本研究条件的急性、重复、双目标细胞组合。 |
| GSE276202 | 1月stroke；9只鼠，3个区域文库 | 9只鼠不是9个独立表达文库，且非急性Sham比较。 |
| GSE171393 | 2 Sham与2 2VO hippocampus文库 | 模型、区域/时窗不同；并不增加病例重复。 |
| GSE254550 | Sham/tMCAO为skull及meninges，brain仅UVB治疗项 | 不能跨脑边界组织拼成病例对照。 |
| GSE279462 | 文章引用过此编号；本次官方SOFT为404 | 仅被引用不能证明当前可下载、可追溯或可纳入。 |
| GSE261494（既有） | D14永久远端MCAO，young Sham 2、stroke 4；aged同样2/4，含两性 | 可以做晚期延伸，但改变急性问题，且每年龄Sham仍2。 |

每个GEO对应的标题、GSM、特征和文件URL均保存在JSON；上述排除并不否定原研究针对其自身问题的价值。

## 新题目与已有方法学工作

不能把“CCC跨样本复现性很少有人研究”当作已证实新颖性。以下4篇是直接需要区分的原始方法/benchmark研究：

1. **Dimitrov 等，2022。**系统拆分7方法与16资源，已表明资源选择、评分/优先级策略及细胞注释影响预测，包含多个数据集及噪声稳健性检查。本稿可区别的对象是同一病理对比的研究间效应转移，以及表达变化与共识排序变化的关系。[原文](https://www.nature.com/articles/s41467-022-30755-0)。
2. **Dimitrov 等，2024，LIANA+。**已有多样本条件比较、不同表达/特异性分数及factorization；用具有多个样本的图谱评价条件分类，并构造随机基线。不能把“每文库运行LIANA”或“多方法共识”本身当创新。[原文](https://www.nature.com/articles/s41556-024-01469-w)。
3. **Cesaro 等，2025，scSeqCommDiff。**差异通讯benchmark已比较CellChat、scDiffCom、MultiNicheNet与LIANA+，使用心肌梗死及多发性硬化的单细胞和匹配空间数据。其重点是差异互作与空间证据，不是本稿的缺血脑两研究效应复现率。[原文](https://doi.org/10.1093/nargab/lqaf084)。
4. **Ku 等，2026。**9种空间CCC工具、9个真实数据集/3项研究，明确量化跨样本和切片Jaccard复现率，并比较原生与共用资源。跨样本通讯复现性已是其正式评价维度；本稿应把增量写成表达效应与排序效应的断裂、整体结构保留null及候选宇宙的控制。[原文](https://doi.org/10.1186/s13059-026-04063-5)。

建议定位为**脑缺血公开队列中的可复现性案例研究与评分层级比较**。若要把标题中的“does not transfer”上升为一般规律，需要更多独立病种/队列与统一候选集合证据；当前可用“limited transfer”或在标题明确是mouse cerebral ischemia范围。对于为什么排序下降，先实证拆分候选资格、背景细胞组成、资源rank母集和评分归一化，再讨论相对优先级与绝对表达量不是同一估计对象；这比笼统宣告某个工具不可复现更有信息。
