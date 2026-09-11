# 论文与可复现分析包

本次交付完成了基于真实公开数据的中文研究稿、5张主图、2张补充图、完整结果表及分析代码。数据检索主要截止于2026年9月9日，交付核验完成于2026年9月10日。

## 先读什么

- `outputs/manuscript_zh.docx`：可编辑中文全文，含摘要、方法、结果、讨论、2个正文表、22篇已核对 DOI 的参考文献及7张图。
- `outputs/manuscript_zh.md`：同一内容的文本版。
- `outputs/figures/`：400 dpi PNG 及可编辑 SVG、矢量 PDF。
- `outputs/tables/`：全部正式结果；TSV 为制表符分隔，`.gz` 为压缩文本。
- `outputs/figure_source_data/`：各图对应的数值及样本记录。
- `work/`：方案、代码、固定资源、来源审计、汇总表达、环境版本与校验记录。

## 本文实际完成了什么

纳入4个研究数据集：两个具有可追溯动物文库的单细胞队列、一个分选细胞队列和一个空间转录组数据集。LIANA 提供相互作用注释，SpatialFeatureExperiment 提供条码几何映射；它们不是额外的独立疾病验证队列。本研究没有为了凑“五个数据库”把注释资源当作生物学重复。

主结果是共同可评估候选中的方向一致性（258/323，79.9%；发现候选中154/182，84.6%）及屏障相关转录异质性。候选评分是自定义的样本级表达共可用性，不是实际 CellChat 概率，也没有声称运行 LIANA 多算法整合。全文保留了分选细胞复核不足、相反基因变化、精确置换检验分辨率不足及环境 RNA 风险。

这是一份有真实结果的完整探索性论文初稿，不能直接称为已完成机制验证或保证录用的稿件。规则式细胞识别、小样本外部队列（3对2）、环境 RNA 未校正、空间切片重复不足是实质性限制。2026年9月已经存在接近主题的多队列整合研究，本文没有宣称首次发现胶质—血管通讯或 SPP1 机制。

## 投稿前需要实际作者完成的事项

核对并填写作者、单位、贡献、资金、利益冲突和本单位要求的伦理声明；确认目标期刊的语言、格式及人工智能辅助披露政策；对原始数据来源、注释选择和统计解释承担最终责任。实际作者未提供的事实没有被代填，也没有代替作者提交稿件、联系期刊或声称已公开注册分析方案。

## 复现方法

建议 Python 3.12。进入解压后的项目目录，创建虚拟环境并安装 `work/requirements_analysis.txt`。相关统计软件版本以 `work/software_versions.json` 为准。当前交付包含已计算的 pseudobulk、检测比例、细胞平均表达、分选原始标准化矩阵和派生空间点数据，常规重算不需要重新下载大型单细胞矩阵。

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r work\requirements_analysis.txt
.\.venv\Scripts\python.exe -X utf8 work\run_pipeline.py
```

默认重算 pseudobulk 差异表达、所有候选评分与下游结果，覆盖包内派生结果。每个数据集各处理组合依次计算，使用1个统计工作线程。若仅需由现有表重建汇总和图文，使用 `--use-cache`。

```powershell
.\.venv\Scripts\python.exe -X utf8 work\run_pipeline.py --use-cache
```

若要从作者提交的单细胞计数矩阵重做质控、注释和汇总，使用 `--from-raw`。程序会创建独立的 `fresh_run_时间戳` 目录，复制固定分析资源，下载公开矩阵后重新计算，不依赖包内已有细胞注释。需要网络、额外数GB磁盘和较长运行时间。这里的“raw”指作者提交的 UMI 矩阵，不包括 FASTQ 重比对。分选数据本来就是作者提供的标准化矩阵，不能转称原始计数。

```powershell
.\.venv\Scripts\python.exe -X utf8 work\run_pipeline.py --from-raw
```

图文生成依赖中文字体；当前稿件使用宋体/黑体及 Times New Roman，图形使用 Arial。不同操作系统字体替换可能改变分页，不应据此认定数值发生变化。Word 已在本次 Windows 环境用 LibreOffice 渲染逐页检查。打包不包含虚拟环境、文档渲染器或1GB级单细胞原矩阵；下载地址及哈希已保存。

## 结果表索引

- `tables/GSE174574/` 与 `tables/GSE245386/`：主分析与3种敏感性的全基因差异、标准化表达、候选全表及逐动物评分。
- `tables/cross_cohort/`：全基因复核、全部候选的同向/反向/表达不足记录和汇总。
- `tables/integrated/`：屏障基因效应、面板逐样本得分、候选组分、留一样本敏感性及探索性相关。
- `tables/sorted_validation/`：6对分选样本的全转录组配对及非配对统计结果。
- `tables/spatial/`：5张切片的质控、真实坐标及描述性相关。
- `tables/identity_marker_audit.tsv`：目标及非目标标记检出审计。
- `tables/validation_checks.json`：190项数值、映射及结果一致性检查。该记录不是新增生物学验证。

不同输出压缩文件的 gzip 时间戳可能改变文件哈希；数值复现应对读取后的表比较。所有输出的小数精度由底层结果表保留，正文为四舍五入展示，未把极小 P 值写成零。
