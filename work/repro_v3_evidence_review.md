# 第三轮独立证据核验

核验完成于 2026-09-11T07:32:38.857253+00:00。仅新增本说明与同名JSON；未改主稿或原计算。

## 已知上游修正

GitHub官方API重新取得：PR #261于 **2026-09-09T07:06:26Z** 合并，完整提交 **d4211373692e7b9c10210488ccb1efe06452b097**。本稿可量化已知修正的下游后果，不能自称首次发现，也不能用较晚预印本追溯取得优先权。仅执行聚合函数与完整升级LIANA是不同操作。 [PR记录](https://github.com/scverse/liana/pull/261)；[固定提交](https://github.com/scverse/liana/commit/d4211373692e7b9c10210488ccb1efe06452b097)。

## 文献与数据来源支持范围

**参考3：Giulia Cesaro，2025-06-19，DOI 10.1093/nargab/lqaf084。** [官方全文](https://academic.oup.com/nargab/article/7/2/lqaf084/8169143)。

该文提出scSeqCommDiff/CClens，比较实验条件间通讯并以配套空间数据检查差异候选的空间共定位支持；可支持“既有差异通讯方法及空间支持评价”的背景。 不能据此宣称本文首先比较条件效应或首先使用空间支持；该DOI是方法研究，不应误写为综述，也不能替本文原两队列提供机制验证。

**参考4：Li-Ting Ku，2026-04-08，DOI 10.1186/s13059-026-04063-5。** [官方全文](https://link.springer.com/article/10.1186/s13059-026-04063-5)。

九种空间通讯工具的模拟及真实数据评价已采用共同LR资源，并以Jaccard衡量跨样本/切片一致性；可支持本文需限定为具体疾病案例和实现后果评估的定位。 DLPFC四张切片来自同一供体，不能写为四个独立生物学重复；重叠、文献通路支持与真实通讯机制是不同证据。

**参考8：Ruolin Zhang，2026-09-03，DOI 10.1002/advs.77547。** [官方全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13542395/)。

正文说明同侧纹状体的Sham及再灌注后0h、24h、7d、14d取样；MCAO阻断60分钟；数据获取段对应GSE332910；Table S1中Sham/24h均列3个生物学重复和3个RNA文库。 未找到逐RNA文库的独立鼠ID或混样规模映射，不能把每个文库明确写成单鼠；原文的NSC机制不构成本文星形—内皮通讯的功能验证。

Zhang全文以[官方PMC Cloud XML](https://pmc-oa-opendata.s3.amazonaws.com/PMC13542395.1/PMC13542395.1.xml)核验，MD5=6e023c481fb842370ff36b69fc633db2；[作者Table S1](https://pmc-oa-opendata.s3.amazonaws.com/PMC13542395.1/ADVS-9999-e77547-s001.xlsx)重新读取，SHA256=c20ec5d550af876a0490b763cdf1ea5a315a15b46894588da4fe2c09fd6e92f7。

[GSE332910官方元数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332910)显示9月9日公开、9月10日更新，关联PMID 42693582；六个所选GSM的RNA标题和处理已复核。平台元数据有T7/Tx内部差异，宜保留MGI/DNBSEQ表述，不能静默选定型号。Series中的泛化免疫措辞不能覆盖GSM处理字段与正文的MCAO/R设计。

## Family C：直接核算200配置

每项都只有observed配置199达到其观测上尾值；严格更大者为0，包含自身的精确并列数为1，1e−12容差下仍为1。因而可写“observed allocation was the unique upper-tail maximum among all 200 allocations”，不宜无方向限定地写“most extreme”。

| 选择 | 统计量 | observed | 次高值（去掉observed） | 更大数 | 并列数（含自身） | 单侧P |
| --- | --- | --- | --- | --- | --- | --- |
| primary | same_direction_all_fraction | 0.798761609907 | 0.780730897010 | 0 | 1 | 0.005 |
| primary | same_direction_nonzero_fraction | 0.798761609907 | 0.780730897010 | 0 | 1 | 0.005 |
| primary | spearman_rho | 0.742788297610 | 0.740173828719 | 0 | 1 | 0.005 |
| reference_singlet | same_direction_all_fraction | 0.757668711656 | 0.739273927393 | 0 | 1 | 0.005 |
| reference_singlet | same_direction_nonzero_fraction | 0.757668711656 | 0.739273927393 | 0 | 1 | 0.005 |
| reference_singlet | spearman_rho | 0.703021713206 | 0.688359437048 | 0 | 1 | 0.005 |

定义为 `P = #{g: T_g >= T_observed - 10^-12}/200`，穷举全部允许配置并包括真实标签，故P=1/200=0.005；不是199次随机模拟，不加Monte Carlo伪计数。若另按距0.5的同向率偏离或距0的rho绝对值衡量，两端共有2个同等极端配置，不能把已报告单侧P改称双侧P。

六个P相同本身不证明高相关。直接检查发现全部400行均无零效应候选，因此每种选择的全分母和非零分母同向率在200配置中完全重合；其余统计向量的Spearman相关为0.961648–0.981791。这证明本数据中的统计依赖，不意味着200个生物学重复或六次独立确认。

## Holm与摘要边界

原六项家族保持不变，Holm校正六项均为0.030。Holm在各边际P值有效时，可在任意依赖下强控制FWER；依赖可能使它保守，但不能把校正值称为实际联合错误率，也不能据此断言真实证据“比0.030更强”。[R官方说明](https://search.r-project.org/R/refmans/stats/html/p.adjust.html)。本次单侧P的有效性仍依赖研究内标签可交换性；事后探索家族也不升级为确认性证据。

Family A：24/30项 nominal P<0.05，但Holm后0项，最小0.15；Family B：12/24项 nominal P<0.05，但Holm后0项，最小0.12。A/B应写未达到各自家族校正显著性；C仅支持其探索性联合标签零基线，不能替代方法间显著差异或外部队列单独复制。

**建议方法措辞：**

We retained the six-entry exploratory family specified in the revision-stage plan and applied Holm correction, which controls the family-wise error rate under arbitrary dependence when the marginal P values are valid. The identical P values reflected the same attainable tail count and did not by themselves establish dependence. Direct inspection showed two pairs of identical concordance statistics and strongly associated allocation vectors, so the six entries were not six independent confirmations; this did not justify redefining the family after observing the results. This plan followed inspection of earlier manuscript results and was not preregistered.

## 四句英文摘要建议

Cross-cohort agreement of ligand–receptor expression effects may differ from agreement of inferred candidate priorities, complicating the selection of reproducible interactions.
We evaluated astrocyte–endothelial candidates in two mouse cerebral ischemia cohorts comprising 11 animal libraries using matched candidate sets, exhaustive animal-label allocations and controlled aggregation analyses, including an independently reproduced, previously merged upstream correction.
On 187 common candidates, expression effects agreed for 80.2% versus 56.1% for native LIANA 1.10.0 priorities, and ranking each distinct score column once increased priority correlation to 0.519; dynamic-gate expression concordance was 79.9% versus a 50.0% joint-null median (one-sided exact P=0.005; six-test Holm P=0.030, exploratory family C), whereas fixed-metric and expression-minus-priority families A and B did not reach Holm-adjusted significance (minimum adjusted P=0.15 and 0.12).
These descriptive and exploratory findings do not establish significant between-method superiority and support separately auditing expression effects, software implementations, eligible candidate sets and directional denominators before selecting candidates for validation.

上述第三句保留预期关键数字及A/B、C边界；它是四句压缩版，不引入新检验。若主稿需纳入第三队列和去污染完整结果，应另按期刊字数要求调整，而不能据压缩摘要隐去正文失败边界。

完整来源链接、远程响应校验和、逐配置极端/并列计数、相关矩阵与原文件SHA见 `repro_v3_evidence_review.json`。
