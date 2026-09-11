# 复现附录与补充分析

Sihuan Zhu

本补充材料提供主稿的完整操作方法、数值溯源及辅助生物学分析。各置换家族及诊断输出保持独立标识。

## 复现附录

### 数据和细胞选择

GSE174574纳入3个假手术与3个MCAO后24小时文库；GSE245386纳入新生成的野生型文库，假手术3个、缺血2个，排除Lrg1敲除和重复使用的GSE174574数据。后者为60分钟缺血后24小时再灌注，与自缺血开始计24小时存在约1小时差别。样本身份由GSM、标题、文件名及原论文共同核对；GSE245386的2025年更正纳入来源记录。[1](https://doi.org/10.1177/0271678X211026770)[2](https://doi.org/10.1186/s12974-023-02941-4)[3](https://doi.org/10.1186/s12974-025-03610-4)

更新检索取得145条GEO GDS命中，对21个新增候选提取完整SOFT并核对原论文；检索式、全部命中及排除理由存档。GSE332910提供同侧纹状体单核数据，补充表S1明确报告Sham和24小时组各3个生物学重复、3个RNA文库。作者采用60分钟MCAO及24小时再灌注；每文库单鼠或混样规模未详述。该队列作为新增描述性扩展保留，不与原3/3和2/3动物配置混入同一置换分布。[4](https://doi.org/10.1002/advs.77547) 该更新检索不作为穷尽性系统综述。

作者提交矩阵按文库质控：基因数200–6,000、UMI≥500、线粒体比例≤20%，重名基因符号合并。固定谱系标记规则区分星形胶质细胞、内皮细胞和其他大谱系，保留未分配细胞。两个原队列另用独立训练的CellTypist全鼠脑模型审计全部质控细胞，并以皮层/海马模型作参考选择检查；逐文库Scrublet标记双细胞，种子20260910。[5](https://doi.org/10.1126/science.abl5197)[6](https://doi.org/10.1016/j.cell.2021.04.021)[7](https://doi.org/10.1038/s41586-023-06812-z)[8](https://doi.org/10.1016/j.cels.2018.11.005) 主要对比包括原规则选择与“原规则和全脑参考一致、参考得分≥0.5且Scrublet未标记”的选择。参考标签未按通讯结果调整，模型得分不解释为校准后的正确概率。完整标记面板、竞争类别、逐细胞分数及版本保存在归档配置与标签表。

规则标签在每文库内计算：标记基因log1p(CP10k)按所有质控细胞标准化（总体标准差下限0.25），截断至−3到10，再对每谱系面板取平均。最高平均分≥0.5、与第二名差≥0.25且该谱系至少2个标记有原始计数时分配标签，否则保留未分配。星形面板为Aldh1l1、Slc1a3、Slc1a2、Gja1、Glul、Sox9、Aldoc；内皮面板为Pecam1、Cdh5、Esam、Tek、Erg、Vwf、Kdr。其余13个竞争谱系的固定面板同列于marker_panel.json；这些标签是规则操作定义。

### 动物级评分和共同候选

固定使用mouseconsensus的3,989条资源记录，在星形胶质细胞→内皮细胞及相反方向计算每个动物的候选分数。表达共可用性为配体与受体细胞平均log1p(CP10k)之积的平方根，复合物取亚基均值中的最小值。原组别资格规则要求全部亚基在至少一组的至少2个动物中检出比例≥10%，发送端与接收端分别判断；5%及20%结果保留作敏感性。效应为MCAO动物平均分数减去Sham动物平均分数。

实际LIANA 1.10.0通过rank_aggregate逐动物运行CellPhoneDB、Connectome、log2FC、NATMI及SingleCellSignalR，保持相同资源，expr_prop=0.1、min_cells=30、n_perms=100、seed=20260910、n_jobs=1、return_all_lrs=False。[9](https://doi.org/10.1038/s41467-022-30755-0)[10](https://doi.org/10.1038/s41556-024-01469-w) 计算每文库所有已分配且≥30细胞的大谱系之间的完整网络，再提取目标两方向。参考支持版本对所有谱系去除双细胞，目标谱系另加参考一致筛选。输入按每细胞全部基因总量归一化。

评分按解释分为表达取向与相对排序两族。共可用性、lr_means、expr_prod与lrscore是相关的表达汇总；Connectome和NATMI共享expr_prod，不构成两份独立验证。五个方法接口在幅度共识中贡献lr_means、expr_prod与lrscore三个独特分数。magnitude_rank是相对名次经稳健秩聚合得到的分数，取值上限为1，并非简单百分位名次；将其转换为1减该分数，使较大值表示较高相对优先级。LIANA的细胞标签置换P值描述文库内背景，不替代动物级疾病检验。缺失条目代表未获得合格分数，始终不以零填补。

主比较使用两研究所有11个动物均完整观测、且所有比较指标有限的相同有向候选集合。报告同向数除以全部共同候选数，以及除以两队列效应均非零的候选数，两种分母并列；绝对效应≤10^-12按数值零处理。Spearman相关覆盖全部共同候选。另按同一共同集合对每动物表达分数作平均并列秩的百分位转换，再计算疾病效应。这一固定集合重排用于观察表示变换，不假装重建了完整LIANA网络。

在GSE332910的6个文库上，按原primary标记规则保留通过QC且每类至少30个核的背景细胞，逐文库运行LIANA 1.10.0完整细胞类型网络。归一化使用全部提交基因的原始UMI总量，采用与原队列相同的float32 log1p(CP10k)、mouseconsensus资源、10%检出阈值、100次内部置换及固定种子。所有归一化稀疏矩阵的存储值均为正；每类pseudobulk与独立分块聚合一致。每个原生magnitude rank先从完整网络精确重建，随后在同一完整网络内，对每个不同的magnitude分数列仅进行一次排名，并对这些排名实施RRA聚合，得到明确标注的unique-column实现诊断。比较使用在全部17个文库、全部相关LIANA分数上均完整的相同99条目标候选，并在这一集合内计算custom、lr_means、expr_prod、lrscore、原生及诊断优先级的缺血减Sham效应；排序始终在各文库完整网络上完成后再取目标子集。另计算原两队列在相同99条候选上的参照。数值TSV以round_trip精度读取，全部比较为描述性；第三队列仅使用raw/primary版本，未作参考支持或去污染版本扩展。

### 保留依赖结构的标签零基线

为保留共享配体、受体及候选间依赖，将每个动物的完整多变量分数向量固定，分别在两研究内穷举3/3和2/3处理标签，共20×10=200种联合配置；不置换基因、候选边或单细胞。由于归一化、参考选择及文库内LIANA评分未用疾病标签，重标记时可复用其逐动物分数。对每种配置计算两种同向率及Spearman相关。原10%组别门槛另在每次置换中重新筛选，因此该基线也包含候选选择步骤。

单侧精确P为200个配置中统计量不小于观测值的比例，含原配置，不增加蒙特卡洛伪重复。报告零分布中位数及2.5%–97.5%分位范围；后者是随机标签参照区间，不是观测效应置信区间。计算前的修订计划规定三组Holm校正：固定集合的2种选择×5指标×3统计共30项、表达减共识优先级的24项差值统计，以及动态门槛的6项探索统计。全部结果均保留。标签无关联零假设下的差值检验不等于检验算法具有相同可复现性；该修订是在已知旧稿结果后制定，未作前瞻注册。

在联合零基线完成后，另作条件置换审计：固定一个队列的观测标签，仅穷举另一队列的20或10种配置，保持相同评分及各自资格规则。该事后诊断用于区分“两队列均无标签关联”的联合零假设与单个队列的可分辨证据；不改变原Holm家族，不将200种联合配置解释为200个独立动物或外部队列样本量的增加。

### 聚合实现和排名母集诊断

实现审计发现所用版本按方法遍历分数列并原地重排，Connectome与NATMI共用的expr_prod列会经过两次排序。用4行确定性输入直接调用已安装的官方函数验证这一行为，并保存对应上游标签源码和哈希。随后逐动物重算完整网络：先验证目标分数及原生共识与既存结果一致，再在同一完整网络上将每个独特的幅度分数列只排序一次，调用相同的稳健秩聚合函数生成诊断共识。[11](https://doi.org/10.1093/bioinformatics/btr709) 不修改安装包；所有输出明确区分原生版本与独特分数列诊断，后者不声称为官方修复版本或生物学真值。该追加比较为描述性实现敏感性，不新增确认性P值。

进一步固定所有11个文库均具有三个完整独特表达分数的全网络边交集，在该全局母集内重新聚合独特分数列，再提取原固定目标边。此步骤只干预排名母集，沿用已计算的文库内表达分数，并未固定上游细胞组成或表达背景。记录每文库边数、细胞类型数、聚合饱和、零效应及名次冗余。诊断、零效应阈值和源代码均存档。

另从[LIANA PR #261](https://github.com/scverse/liana/pull/261)的固定合并提交d4211373692e7b9c10210488ccb1efe06452b097取得聚合函数，保存源码和SHA256，并在22份原始完整网络上独立执行，以逐值比较独特分数列诊断。仅核验聚合函数，未安装或重跑完整新版流程；原资源、显式expr_prop=0.1和表达分数保持不变，避免把同一上游变更中的默认检测门槛变化混入比较。

### 环境RNA与数值收敛敏感性

将 DecontX 作为计算性环境 RNA 敏感性分析逐库运行（R 4.6.1、Bioconductor 3.23、decontX 1.10.0）。输入为 GSE174574 和 GSE245386 的 11 个文库中全部原始质控合格细胞的整数 UMI 计数；重复基因符号在模型拟合前合并求和。固定使用独立全脑参考注释的 whole_brain_broad 标签作为 z，不按目标候选或拟合结果重新选择细胞。参数为 background=NULL、batch=NULL、seed=20260911、delta=c(10,10)、estimateDelta=TRUE、convergence=0.001、iterLogLik=10，数值线程设为 1。初始 maxIter=500；保留全部初始输出，并按统一数值规则仅对未达到收敛阈值的文库用相同输入、分群和随机种子将 maxIter 延长至 2000 重新拟合。最终来源仅依据数值收敛选择：初始收敛者使用默认输出，未收敛者使用延长拟合输出；选择不依据候选通讯结果。校正矩阵保存为 float64 的分数估计，提供 cells×genes CSR NPZ 与 genes×cells MatrixMarket；不取整、不将其作为原始 UMI 或输入 DESeq2，并保留基因及细胞顺序、校正后每细胞总量和来源 SHA256。该模型在缺少空液滴实测背景时从已过滤细胞推断非本细胞来源成分，依赖所给细胞群标签；参考身份不确定及细胞群间共享的真实表达可能影响分离结果。污染比例因此属于模型估计而非实测污染率；该分析不直接测量配体释放、受体活化或屏障功能。 [12](https://doi.org/10.1186/s13059-020-1950-6)

连续非负校正计数保留原精度及细胞/基因顺序，在固定细胞选择下重新计算表达分数和完整网络。归一化分母为校正后每细胞全部基因总和。LIANA沿用float32归一化表示，记录精度转换前后正值数，并在计算稀疏检出比例前移除显式数值零，以实施>0检测规则；自定义表达评分另输出估计计数≥1的资格敏感性，以显示连续估计与原始整数检测的区别。零总量细胞触发审计，不静默套用原UMI分母。校正估计不当作原始UMI输入DESeq2，原始pseudobulk及其全部检验保留存档。

### 分选和空间补充分析

GSE163752使用6对独立样本池，每池3–4只鼠，对照为对侧半球。仅采用元数据映射一致的星形胶质与内皮样本，根据提交的重复编号配对，对作者标准化表达作log2(x+1)，报告缺血侧减对侧的平均差及点式95% t置信区间（自由度5）。[13](https://doi.org/10.1007/s00401-022-02452-1) 不以这些区间作同时显著性判定，不把对侧等同Sham；原全转录组配对/Welch及BH结果完整保留，展示方式变化不更改先前检验结果。

屏障相关转录效应作为细胞状态背景，原始计数按动物–细胞类型汇总后进行PyDESeq2分析，并报告固定面板全部基因。[14](https://doi.org/10.1186/s13059-014-0550-8)[15](https://doi.org/10.1093/bioinformatics/btad547)[16](https://doi.org/10.1038/nature25739)[17](https://doi.org/10.1038/s41593-019-0497-x) GSE233814的5个切片记录仅保留一张真实坐标下的补充表达图；不将空间点当动物，不分配未经确认的梗死区域，也不以描述性相关支持中心结论。[18](https://doi.org/10.1073/pnas.2404203121)

### 完整拟合与数值诊断

完成了 11 个文库共 115,355 个原始质控合格细胞的 DecontX 拟合。10 个文库在默认 500 次迭代上限内达到 0.001 收敛阈值；唯一未收敛文库 GSM5319992 在相同输入和设置、仅将上限改为 2000 后于第 637 次迭代达到阈值（末次最大参数变化 0.0009979）。最终采用 10 份默认结果及该文库的延长拟合结果，11 个最终来源均达到预设数值阈值，全部默认结果仍保留。合并细胞后的估计污染比例中位数为 6.13%，各文库中位数范围为 3.98%–8.16%；合并统计仅作描述，不把细胞作为独立生物学重复。校正后每细胞总量均大于零，11 个最终矩阵均为 float64 分数估计。GSM5319992 延长前后估计污染中位数由 6.32% 变为 6.14%，计数矩阵绝对差总和相当于其原始总 UMI 的 0.3565%；不过单细胞污染估计最大绝对差为 37.54 个百分点，故总体差异较小不能替代逐细胞或下游结果核查。这些结果描述模型拟合与数值敏感性，不单独支持特定通讯候选或屏障调控机制。

### 独立CellChat方法与溯源

该探索性扩展在R 4.6.1中使用完整原生 CellChat R 包2.2.0.9001，固定官方上游提交为75253cd0c9e68410e6e721a6d3a0419a1d7e358f。Windows二进制包来自 r-universe 对官方源码的第三方构建；构建来源、已安装核心函数及完整原生小鼠数据库均与固定源码进行了核对。原11库分别分析：GSE174574为3个Sham与3个MCAO文库，GSE245386为3个Sham与2个MCAO文库。采用固定的主要QC细胞标签，并保留每库至少30个细胞的全部背景细胞类型。本扩展未纳入第三队列、参考筛选或DecontX校正后的输入。[19](https://doi.org/10.1038/s41596-024-01045-4)[20](https://doi.org/10.1038/s41467-021-21246-9)

展开完整复合物亚基并按规范顺序排列，未作跨物种转换或别名替换。在原生小鼠蛋白信号目录与mouseconsensus共享的配体—受体定义中，1,548项具有唯一的原生行映射，且其配体、受体和所引用辅助因子基因在全部原始矩阵中均有覆盖。保留原生复合物和辅助因子定义。每细胞按log1p(10,000×计数/全部基因UMI总量)标准化，复现既存LIANA的float32运算后，无损提升至float64作为R输入。仅在完成标准化后导出信号基因，未使用该基因子集的计数总量作为分母。

每库创建原生CellChat对象并指定控制目录，调用subsetData后，以显式LR.use、type="triMean"、raw.use=TRUE、population.size=FALSE、nboot=100、seed.use=20260911、Kh=0.5和n=1调用computeCommunProb。逐库顺序执行，不作PPI投影或空间距离建模。显式LR.use绕过了每库过表达互作预筛选，因此属于控制配置，非默认的完整分析流程。原生函数还按每库信号表达最大值进行缩放。包中的net$prob虽命名为probability，实为推断通讯强度，并非校准后的生物学概率。保留全部零得分及库内细胞标签检验P值；这些P值既不用于筛选本次比较，也不作为动物层面的疾病检验。

在读取新增跨队列CellChat结果前，比较网络固定为控制资源中、在所有11库均有完整LIANA分量得分的1,222个来源细胞—靶细胞—互作条目；这些条目也均存在于各原生输出。该选择使比较以LIANA检测资格为条件。每库由研究者计算的强度百分位为（原生强度的升序平均秩−1）/(1,222−1)，数值越大表示相对强度越高，并非CellChat原生共识排名。在完全相同的1,222条上重算LIANA单列RRA，lr_means、expr_prod及lrscore各使用一次，并恢复既存分量的float32数值表示。RRA优先级定义为1−RRA，未复用其他网络背景的既存排名。相同32个目标条目也按原公式计算表达共可用性。

各指标的队列对比均为文库层面MCAO均值减Sham均值。分别报告Spearman相关系数、以全部32条为分母的同向计数、以两队列效应均非零条目为分母的同向计数，以及零效应数。应用并审计原有10^−12绝对效应容差。保留原生RDS，通过17位有效数字的无损导出与原RDS核对后，再以可往返精度读取。原A/B/C置换检验P值不适用于改变后的指标和候选集合，本扩展未进行新的疾病效应推断检验。方向与辅助因子分层仅作描述。单独的完整原生资源敏感性分析未运行，因此本分析不评估CellChat在其完整原生目录中的性能。

### 来源与发布结构

仓库路径与哈希列于FILE_INDEX.tsv，图源数据由FIGURE_SOURCE_INDEX.tsv映射。ENVIRONMENT.md区分实际验证环境和恢复方案；RUNNING.md区分缓存重建、原始推断与固定提交聚合函数重放。公开版本标识列于主稿代码获取段。

## 受控独立框架结果

资源审计保留原187个有向候选中的32条（17.1%）：星形胶质细胞→内皮细胞8条、内皮细胞→星形胶质细胞24条，共涉及29个不同的配体—受体定义。其余155条在原生CellChatDB中缺少完整亚基对的精确匹配；在共同全网络取交集时未再丢失目标条目。这些排除源于资源定义，并非CellChat预测失败。保留条目中31条不含原生辅助因子，1条含辅助因子，无法据此计算有意义的辅助因子分层相关性。

在共同32条中，表达共可用性同向23/32（71.9%；ρ=0.565），原生强度同向10/32（31.3%；ρ=0.549），研究者计算的百分位同向16/32（50.0%；ρ=0.364），相同网络的LIANA RRA优先级同向13/32（40.6%；ρ=0.451）。原生强度在发现与外部队列分别有19与3条零效应，两队列均非零的13条中同向10/13（76.9%）。LIANA RRA分别有10与12条零效应，两队列均非零的19条中同向13/19（68.4%）。表达共可用性和百分位均无零队列效应。对这四个指标，两队列均无非零对比因10^−12容差而被置零。

星形胶质细胞→内皮细胞Lamb2–Dag1、Lamb2–Itga6_Itgb1以及内皮细胞→星形胶质细胞Entpd1–Adora2b，在全部11库中的原生强度均为零，但各自的百分位对比在发现队列约为+0.014060、外部队列约为−0.024775。每库1,222个网络条目中有160–862个零强度，平均秩赋予零并列块的百分位下限为0.065111–0.352580。因此，这些百分位变化反映相对背景改变，而其原生零强度并未改变。

分方向看，表达共可用性在星形胶质细胞→内皮细胞方向同向4/8（ρ=0.071），反向同向19/24（ρ=0.696）；原生强度分别同向1/8（非零效应中为1/4；ρ=−0.230）和9/24（非零效应中为9/9；ρ=0.817）。这些小且不均衡的子集仅作描述，不建立方向间或方法间差异。完整得分、效应、资源流失及并列秩审计与补充图S6一并归档。

## 补充生物学分析

### 分选效应与空间表达提供补充背景

分选数据以每个样本池的配对表达差展示，而不以全转录组显著条目数判断复现。11个既定示例组分中，星形胶质Spp1的平均log2(x+1)差为3.376（点式95% CI 0.653–6.100），Timp3为−0.686（−1.096至−0.275），内皮Plat为−0.043（−0.809至0.722）；同一基因的所有6对差值均可查看（补充图S4）。全部49,914个细胞类型–基因描述效应及区间已另存，与旧检验集合相应效应/区间的最大绝对误差小于2×10^-15。

原始pseudobulk中Slco1c1下调、Angpt2上调在原两队列及细胞选择敏感性中保持，其他屏障相关基因有异质性，作为细胞状态背景保留（补充图S3）。空间表达仅列为单张补充图S5，完整切片描述表存档，不参与表达与排序可复现性之间的中心比较。




## 补充图注

补充图S1　独立参考身份与细胞保留。A–D：各位置为一个动物文库，圆点为全脑标签一致率，菱形为皮层/海马一致率，空心方块为全脑一致、最高sigmoid得分≥0.5且Scrublet未标记的保留率，分母为原规则目标细胞。E、F：先计算每动物原标签行内的全脑参考标签比例，再等权平均，格内为四舍五入百分数。最佳匹配参考没有未分配输出；参考一致率不是测得的身份准确率。

补充图S2　各文库Scrublet观测与模拟分数。统一分箱，纵轴为相应观测或模拟总数归一化后的箱内比例；蓝色实线为观测细胞，橙色虚线为模拟双细胞，黑色点划线为自动阈值。GSM7841720的阈值高于全部观测分数，故标记数为零；未手动调阈值。这是技术质控，无细胞层面的疾病推断。

补充图S3　屏障相关单基因的细胞选择敏感性。各格为相应原队列和选择中的内皮细胞动物级PyDESeq2 MCAO相对Sham log2倍数变化。空心圆为相应内皮转录组BH校正P<0.05；色阶截于±3，完整值存档。白线分隔连接、转运、激活及血管生成面板。同队列四列复用相同动物，属于选择敏感性；此图使用原始UMI pseudobulk，未将连续DecontX估计送入DESeq2。

补充图S4　分选样本的配对效应与区间。A为星形胶质组分，B为内皮组分；彩色点为六个独立样本池各自的缺血侧减对侧log2(作者标准化值+1)差，黑点及线为均值和点式95% t置信区间（自由度5）。每池含3–4只鼠，对侧并非Sham。11个组分来自既定示例，未按本图显著性重新选择；所有49,914项描述结果及此前全部P值另存，不以区间作同时显著性判断。

补充图S5　候选组分空间表达背景。每点为质控合格Visium空间点，上行为Control切片（2,416点），下行为Day 1切片（2,498点），A–D为Spp1、Itga5、Ptn、Ptprz1。坐标来自原始alignment JSON与已核验条码映射；同一基因上下图使用相同色阶，上限为两图合并表达的第99百分位。没有组织学区域标签，未指定核心或半暗带；空间点混合多类细胞。本图仅提供表达背景，5个切片记录的全部描述表保留存档。

补充图S6. 控制资源范围的原生CellChat比较。各面板分别展示队列对比：A，表达共可用性；B，CellChat原生推断强度（net$prob）；C，研究者计算的强度百分位；D，在相同1,222条网络背景中重算的LIANA单列RRA优先级。横轴为GSE174574、纵轴为GSE245386，均为文库层面MCAO均值减Sham均值。各面板保留相同32个有向候选。蓝色圆点表示8个星形胶质细胞→内皮细胞条目，橙色三角表示24个反向条目。空心符号在各面板均标记原生强度在全部11库中为零的相同3条；重叠条目可能位于相同坐标。图内标注Spearman ρ、以全部32条为分母及以两队列均非零效应为分母的同向计数，以及任一队列效应为零的条目数。若限定两队列效应均非零，原生强度同向10/13，LIANA RRA同向13/19。C中的空心点显示原生零强度不变时，网络背景改变仍可使零并列块的百分位位置移动。该百分位由研究者计算，并非CellChat原生共识输出。以上比较均为描述性分析。

## 参考文献

1. Kai Zheng, Lingmin Lin, Wei Jiang, Lin Chen, Xiyue Zhang, Qian Zhang, et al. Single-cell RNA-seq reveals the transcriptional landscape in ischemic stroke. Journal of Cerebral Blood Flow & Metabolism. 2022;42:56-73. [doi:10.1177/0271678X211026770](https://doi.org/10.1177/0271678X211026770)

2. Zhaohui Ruan, Guosheng Cao, Yisong Qian, Longsheng Fu, Jinfang Hu, Tiantian Xu, et al. Single-cell RNA sequencing unveils Lrg1's role in cerebral ischemia‒reperfusion injury by modulating various cells. Journal of Neuroinflammation. 2023;20:285. [doi:10.1186/s12974-023-02941-4](https://doi.org/10.1186/s12974-023-02941-4)

3. Zhaohui Ruan, Guosheng Cao, Yisong Qian, Longsheng Fu, Jinfang Hu, Tiantian Xu, et al. Correction: Single-cell RNA sequencing unveils Lrg1’s role in cerebral ischemia‒reperfusion injury by modulating various cells. Journal of Neuroinflammation. 2025;22:269. [doi:10.1186/s12974-025-03610-4](https://doi.org/10.1186/s12974-025-03610-4)

4. Ruolin Zhang, Chang Liu, Yuneng Zhou, Kaichen Zhao, Muyang Li, Bingcheng Cai, et al. Nrsn1–Smarcc1 Coupling Regulates Neural Stem Cell Differentiation and Chronic‐phase Recovery After Ischemic Stroke. Advanced Science. 2026;e77547. [doi:10.1002/advs.77547](https://doi.org/10.1002/advs.77547)

5. C. Domínguez Conde, C. Xu, L. B. Jarvis, D. B. Rainbow, S. B. Wells, T. Gomes, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science. 2022;376:eabl5197. [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)

6. Zizhen Yao, Cindy T.J. van Velthoven, Thuc Nghi Nguyen, Jeff Goldy, Adriana E. Sedeno-Cortes, Fahimeh Baftizadeh, et al. A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation. Cell. 2021;184:3222-3241.e26. [doi:10.1016/j.cell.2021.04.021](https://doi.org/10.1016/j.cell.2021.04.021)

7. Zizhen Yao, Cindy T. J. van Velthoven, Michael Kunst, Meng Zhang, Delissa McMillen, Changkyu Lee, et al. A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. Nature. 2023;624:317-332. [doi:10.1038/s41586-023-06812-z](https://doi.org/10.1038/s41586-023-06812-z)

8. Samuel L. Wolock, Romain Lopez, Allon M. Klein. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. Cell Systems. 2019;8:281-291.e9. [doi:10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)

9. Daniel Dimitrov, Dénes Türei, Martin Garrido-Rodriguez, Paul L. Burmedi, James S. Nagai, Charlotte Boys, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nature Communications. 2022;13:3224. [doi:10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0)

10. Daniel Dimitrov, Philipp Sven Lars Schäfer, Elias Farr, Pablo Rodriguez-Mier, Sebastian Lobentanzer, Pau Badia-i-Mompel, et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. Nature Cell Biology. 2024;26:1613-1622. [doi:10.1038/s41556-024-01469-w](https://doi.org/10.1038/s41556-024-01469-w)

11. Raivo Kolde, Sven Laur, Priit Adler, Jaak Vilo. Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics. 2012;28:573-580. [doi:10.1093/bioinformatics/btr709](https://doi.org/10.1093/bioinformatics/btr709)

12. Shiyi Yang, Sean E. Corbett, Yusuke Koga, Zhe Wang, W Evan Johnson, Masanao Yajima, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biology. 2020;21:57. [doi:10.1186/s13059-020-1950-6](https://doi.org/10.1186/s13059-020-1950-6)

13. Daniel Spitzer, Sylvaine Guérit, Tim Puetz, Maryam I. Khel, Moritz Armbrust, Maika Dunst, et al. Profiling the neurovascular unit unveils detrimental effects of osteopontin on the blood–brain barrier in acute ischemic stroke. Acta Neuropathologica. 2022;144:305-337. [doi:10.1007/s00401-022-02452-1](https://doi.org/10.1007/s00401-022-02452-1)

14. Michael I Love, Wolfgang Huber, Simon Anders. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15:550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

15. Boris Muzellec, Maria Teleńczuk, Vincent Cabeli, Mathieu Andreux. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics. 2023;39:btad547. [doi:10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547)

16. Michael Vanlandewijck, Liqun He, Maarja Andaloussi Mäe, Johanna Andrae, Koji Ando, Francesca Del Gaudio, et al. A molecular atlas of cell types and zonation in the brain vasculature. Nature. 2018;554:475-480. [doi:10.1038/nature25739](https://doi.org/10.1038/nature25739)

17. Roeben Nocon Munji, Allison Luen Soung, Geoffrey Aaron Weiner, Fabien Sohet, Bridgette Deanne Semple, Alpa Trivedi, et al. Profiling the mouse brain endothelial transcriptome in health and disease models reveals a core blood–brain barrier dysfunction module. Nature Neuroscience. 2019;22:1892-1902. [doi:10.1038/s41593-019-0497-x](https://doi.org/10.1038/s41593-019-0497-x)

18. Daniel Zucha, Pavel Abaffy, Denisa Kirdajova, Daniel Jirak, Mikael Kubista, Miroslava Anderova, et al. Spatiotemporal transcriptomic map of glial cell response in a mouse model of acute brain ischemia. Proceedings of the National Academy of Sciences. 2024;121:e2404203121. [doi:10.1073/pnas.2404203121](https://doi.org/10.1073/pnas.2404203121)

19. Suoqin Jin, Maksim V. Plikus, Qing Nie. CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nature Protocols. 2025;20:180–219. [doi:10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)

20. Suoqin Jin, Christian F. Guerrero-Juarez, Lihua Zhang, Ivan Chang, Raul Ramos, Chen-Hsiang Kuan, et al. Inference and analysis of cell-cell communication using CellChat. Nature Communications. 2021;12:1088. [doi:10.1038/s41467-021-21246-9](https://doi.org/10.1038/s41467-021-21246-9)
