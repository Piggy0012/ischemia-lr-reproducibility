# 独立细胞注释、双细胞及污染敏感性修订建议

核对日期：2026-09-10。只读了现有 `stream_annotate.py` 和中文稿件；未改动主分析代码、标签或结果。本文件是实施前方法建议，不是新增分析结果。

## 结论与优先级

优先实施 **CellTypist 成年全鼠脑固定外部模型的分批预测 + 逐动物文库 Scrublet + 分组呈现注释/污染不确定性**。这三个步骤在 Windows Python 和约 8 GB 内存条件下具有可行性。建议先运行一个对照和一个缺血文库试验并记录峰值内存；禁止一次合并全部矩阵。模型由独立健康图谱训练，不能将其标签当真值，更不能用本研究规则标签再训练并称为独立验证。

CellTypist 建模原理、输入和预测见[官方代码](https://github.com/Teichlab/celltypist)及[官方预测说明](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html)。本任务新增的参数属于事后修订，应在修订记录中冻结后再看下游结果。

## 已实际取得的外部模型

模型从 [CellTypist 官方目录](https://celltypist.cog.sanger.ac.uk/models/models.json) 下载，保存在 `work/identity_reference_review/`。完整来源、版本、文件 SHA-256 见 `manifest.json`；全部标签、参考基因及本地匹配审计见 `model_inspection.json`。不需要下载完整百万细胞图谱。

| 模型 | 实际下载体积 | 类别/特征基因 | GSE174574 匹配 | GSE245386 匹配 | 建议角色 |
| --- | ---: | ---: | ---: | ---: | --- |
| Mouse_Whole_Brain.pkl v1 | 7,672,163 bytes | 334 / 5,596 | 5,070（90.60%） | 5,594（99.96%） | 主要外部参考 |
| Mouse_Isocortex_Hippocampus.pkl v1 | 680,650 bytes | 42 / 3,383 | 3,268（96.60%） | 3,327（98.34%） | 参考选择敏感性 |

以上数字来自本地下载文件和11份作者 feature 表的实际检查，而非估算；匹配指模型基因存在于提交 feature 表，不代表在每个细胞中实际检出。两个模型分别对应 Yao 等的[成年全鼠脑图谱](https://doi.org/10.1038/s41586-023-06812-z)与[成年皮层/海马图谱](https://doi.org/10.1016/j.cell.2021.04.021)。两者来自同一研究体系且可能共享图谱构建背景，不应称作两个完全独立的验证实验。

全脑模型含 `317 Astro-CB NN`、`318 Astro-NT NN`、`319 Astro-TE NN`、`320 Astro-OLF NN`、`333 Endo NN`，以及 Microglia、BAM、Monocytes、DC、Lymphoid。模型另有 Astroependymal、Bergmann、Tanycyte 等标签，不应不加判断一并归入 Astro。42 类模型含 Astro、Endo、Micro-PVM、SMC-Peri、VLMC，但没有独立单核、淋巴或 OPC 标签；因此不适合作为缺血浸润背景下唯一的分类参考。

## 实施细节及损伤状态偏差

1. 对每个文库所有质控合格细胞预测，避免只预测原规则已选中的 Astro/Endo 导致分母选择偏差。使用作者完整计数计算每细胞总 UMI，重复符号先合并，然后 log1p(CP10k)；不可只用模型基因重新计算归一化分母，不作人体同源转换或全转大写。按1,000–2,000个细胞分批；全脑模型5,596特征的2,000行 float64 矩阵约89.5 MB，另需计入原始稀疏矩阵和中间副本。
2. 保留完整细粒度最佳标签、最高 decision score、sigmoid confidence、不同大谱系的次高分、未分配和多标签状态。第一版不使用 majority voting，避免聚类传播强制覆盖少数损伤细胞。预测使用官方 API；如实现内存有界的固定系数计算，必须与官方 API 在同一输入上逐值比较后方可采用。
3. 主输出为模型最佳标签的描述性混淆表；另预先指定较严格集合，例如 `p_top >= 0.5` 且不同大谱系的概率差 >=0.1。该阈值只是分析规则，不能叫经本任务校准的90%准确率。CellTypist 的类别 sigmoid 分值不保证总和为1，合并 Astro 子类时不能直接把概率相加。严格阈值下未确认者保留为不确定。
4. 分别呈现每只动物及处理组的原规则/外部模型一致率、不确定率、各 Astro 区域标签构成和 UMI/线粒体分布。参考模型从健康组织学习，反应性细胞可能降低置信度；若缺血组损失更多细胞，不能解释成“缺血组更脏”而直接删除。
5. 主要下游敏感性至少有两条：(a) 仅使用外部模型定义的目标细胞重新聚合；(b) 原规则与模型一致且高置信度的交集重新聚合。二者仍按动物、细胞类型 pseudobulk；每样本门槛保持30细胞，细胞不足就报告该比较不可评估，不补零或降低门槛救结果。原主分析应保留，避免事后用最有利的标签替换。
6. 两队列模型基因覆盖不同，应增加固定共同特征集合预测敏感性，并检查 Astro/Endo 最高正权重基因的覆盖。比较两种特征方案的标签及动物效应；不能因总体覆盖90%就断言不存在覆盖偏差。
7. 将模型置信度和 Spp1、髓系转录本等并列审计。Spp1 的细胞来源需依据多标记组合与跨方法一致性，不能以模型给了 Astro 标签视为来源确证。健康参考之间的一致也无法完全解决损伤状态缺失的问题。

当前模型在 sklearn 0.24.1 序列化，而本地为1.9.0；读取出现 `InconsistentVersionWarning`。两个模型的系数、scaler、基因和标签已经成功读取，但尚未完成推断兼容性试验。应保留警告、检查系数维度/有限性、验证官方预测与公式 `clip_upper((x-mean)/sd,10) @ coef.T + intercept` 一致，不能把成功 unpickle 当成运行已验证。公式顺序、特征缺失处理和标签映射必须与[官方 classifier 源码](https://celltypist.readthedocs.io/en/latest/_modules/celltypist/classifier.html)对齐。

## 双细胞分析

[Scanpy 的 Scrublet 接口](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html)接收未归一化计数，按单个文库预测。每只动物文库单独处理完整可用细胞群，不要先只保留 Astro/Endo，否则缺少合成异型双细胞的来源。保留当前低质量过滤，双细胞分析应置于目标身份筛选之前。建议固定随机种子0、30 PCs、模拟比例2、默认预期双细胞率0.05；缺乏装载信息时该值仅是假设，补0.03及0.10敏感性，不声称为测量的真实双细胞率。

保存 observed/simulated 分数直方图、自动阈值、预测比例以及剔除前后目标细胞数。检查自动阈值和模拟分布是否合理；若没有可解释的分隔，应保留连续分数并明确该文库判定不稳，不为得到预期比例随意画阈值。建议同时提供去除 Scrublet 标记细胞后的原规则标签及外部模型标签两类重分析。Scrublet 不检出所有同型双细胞，也不能替代 RNA 污染校正；作者提交矩阵可能已经过筛选，不能称本步骤完整恢复了所有原始液滴中的双细胞。原方法及执行建议见 [Scrublet 官方仓库](https://github.com/AllonKleinLab/scrublet)。

## 环境 RNA：可做敏感性，不能假称已完成校正

本地所谓原始 UMI 是作者提交的细胞矩阵，并未证明含空液滴。标准 SoupX 需要结合背景液滴估计污染谱；此时不能把低质量细胞或所有细胞的平均表达冒称已测空液滴背景。参见 [SoupX 官方仓库](https://github.com/constantAmateur/SoupX)。

现环境可先保留此前强髓系信号排除，叠加独立模型与 Scrublet，并提供每样本目标谱系中髓系标记检测分布、候选表达集中度及剔除后动物级效应。必须称“混合髓系信号/污染敏感性”，不能改写成“已去除环境 RNA”。同一排除规则对 Sham 和 MCAO 一致执行，报告剔除率及可能移除真实炎症反应细胞的偏差。

如果随后允许补 R 工具链，DecontX 是可仅用 filtered cells 运行的成熟方法；其默认从其他细胞群估计污染，对聚类划分和缺失群体敏感，健康/损伤共享反应转录本也可能被过度扣除。必须逐文库运行，分别用无监督宽谱系聚类与外部参考宽标签评估，并将校正结果作为敏感性。[DecontX 官方说明](https://bioconductor.org/packages/release/bioc/vignettes/decontX/inst/doc/decontX.html)明确允许不提供背景矩阵，也明确依赖细胞群标签。没有 R 时，不建议临时重写其贝叶斯算法并沿用 DecontX 名称。

## SingleR 替代路线及为何暂不优先

SingleR 已有[官方 Python 绑定](https://github.com/SingleR-inc/singler-py)，不是只能在 R 使用。官方 celldex 提供 `mouse_rnaseq`、固定版本 `2024-02-26`，含358个分选 bulk RNA-seq参考样本、18大类，包括 Astrocytes、Endothelial cells、Microglia、Monocytes、Granulocytes 等。它适合大类交叉核对，缺少周细胞/OPC等独立类，且 bulk 与单细胞平台存在偏差。[参考说明](https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html)

实际查询 PyPI 显示 singler 0.5.0 当前提供 Linux/macOS wheels，没有 Windows wheel，需C++构建和其依赖；本边界任务未安装也未取得该参考矩阵，故不伪报精确下载体积。可通过 [celldex Python](https://github.com/SingleR-inc/celldex-py) 获取参考，358列的表达矩阵本身预计不构成8GB瓶颈，但Windows安装风险使它不应阻塞现可执行的 CellTypist 方案。也不能将简单参考均值相关自写脚本称作已运行 SingleR。

## 应交付的新增证据

- 模型和软件版本、固定资源哈希、逐文库基因覆盖及各大谱系权重基因覆盖；
- 所有质控细胞的外部参考标签/置信度、原规则与模型的逐动物混淆表；
- 双细胞分数与阈值图、处理组剔除率、样本保留情况；
- 各注释/双细胞/污染敏感性下完整差异表达、全部候选效应和分母流失情况，而非只展示4个示例；
- 修订讨论清楚区分：独立注释支持、参考状态偏差、双细胞敏感性，以及仍未解决的环境 RNA 和细胞来源不确定性。

本次建议无需新增湿实验，也不改变生物学样本量或现有探索性研究定位。
