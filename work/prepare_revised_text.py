"""Preserve the initial manuscript and assemble documented revision insertions."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent
raw=(ROOT/'manuscript_source.md').read_text(encoding='utf-8')
insert=(ROOT/'revision_text_insertions.md').read_text(encoding='utf-8')
def block(start):
    return insert.split(start+'\n\n',1)[1].split('\n## ',1)[0].strip()
raw=raw.replace('分析截止日期：2026年9月9日','初始分析截止日期：2026年9月9日　方法补强完成于2026年9月10日')
raw=raw.replace('### 样本级差异表达',block('## 方法：接在“单细胞质控与细胞身份判定”之后')+'\n\n### 样本级差异表达')
raw=raw.replace('### 屏障相关转录特征',block('## 方法：接在“候选通讯评分及跨队列一致性”之后')+'\n\n### 屏障相关转录特征')
raw=raw.replace('也不是实际执行 LIANA 多算法整合后的秩。','也不是 LIANA 多算法整合后的秩；实际 LIANA 多方法排序另行计算并报告。')
raw=raw.replace('动物级比较报告实际可评估样本和候选分母，并将表达条件变化与评分变化分别记录。','疾病效应比较仅纳入在相应队列全部动物文库中均有合格分数的候选，计算缺血组减假手术组的均值差、Welch 检验及穷举动物标签的双侧精确置换检验；Welch P 值在每个队列、细胞选择版本和评分变量的全部完整观测候选中分别进行 BH 校正。跨队列比较采用两队列共同的完整观测集合；magnitude_rank 转换为1减秩，使正向差值表示相对优先级上升。该完整观测规则及 return_all_lrs=False 在初次接口检查后、查看 LIANA 疾病差异前固定，属于明确记录的技术修订。不同评分变量的效应单位不能直接比较。')
needle='方向一致性本身也需要谨慎理解。'
raw=raw.replace(needle,block('## 讨论：接在首段之后，明确近期创新边界')+'\n\n'+needle)
needle='四个示例具有不同的证据结构。'
raw=raw.replace(needle,block('## 讨论：接在方向一致性及小样本段之后')+'\n\n'+needle)
needle='屏障相关分析说明，'
raw=raw.replace(needle,block('## 讨论：接在 Spp1 段之后')+'\n\n'+needle)
needle='本研究存在若干具体局限。'
raw=raw.replace(needle,block('## 讨论：替换原“尚缺独立标签验证、完整双细胞识别”相关表述')+'\n\n'+needle)
raw=raw.replace('第二，采用透明规则注释可以重现筛选逻辑，但尚缺独立标签验证、完整双细胞识别和环境 RNA 去除；较宽的谱系类别还可能混合内皮区段或反应状态。','第二，已增加外部参考注释及双细胞敏感性，但仍无受损组织的身份真值，也未校正环境 RNA；健康参考的适用范围及自动阈值的样本依赖性限制了该补强的解释，较宽的谱系类别还可能混合内皮区段或反应状态。')
raw=raw.replace('较宽的谱系类别还可能混合内皮区段或反应状态。','较宽的谱系类别还可能混合内皮区段或反应状态。',1)
raw=raw.replace('本研究没有公开注册的分析方案。','新增参考模型、双细胞及 LIANA 的设置和技术调整另见修订方案。本研究没有公开注册的分析方案。')
(ROOT/'manuscript_revised_source.md').write_text(raw,encoding='utf-8')
builder=(ROOT/'build_manuscript.py').read_text(encoding='utf-8')
builder=builder.replace("OUT=ROOT.parent/'outputs';OUT.mkdir(exist_ok=True)","OUT=ROOT.parent/'outputs/revised';OUT.mkdir(parents=True,exist_ok=True)")
builder=builder.replace("ROOT/'manuscript_source.md'","ROOT/'manuscript_revised_source.md'")
(ROOT/'build_revised_manuscript.py').write_text(builder,encoding='utf-8')
print('Prepared revised source and builder; initial manuscript preserved')
