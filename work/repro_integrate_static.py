"""Integrate completed diagnostics without altering underlying results."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
p=ROOT/'repro_main_sections.json'
d=json.loads(p.read_text(encoding='utf-8'))
conditional='在联合零基线完成后，另作条件置换审计：固定一个队列的观测标签，仅穷举另一队列的20或10种配置，保持相同评分及各自资格规则。该事后诊断用于区分“两队列均无标签关联”的联合零假设与单个队列的可分辨证据；不改变原Holm家族，不将200种联合配置解释为200个独立动物或外部队列样本量的增加。'
if conditional not in d['methods_null']:d['methods_null'].append(conditional)
background='进一步固定所有11个文库均具有三个完整独特表达分数的全网络边交集，在该全局母集内重新聚合独特分数列，再提取原固定目标边。此步骤只干预排名母集，沿用已计算的文库内表达分数，并未固定上游细胞组成或表达背景。记录每文库边数、细胞类型数、聚合饱和、零效应及名次冗余。诊断、零效应阈值和源代码均存档。'
if background not in d['methods_implementation']:d['methods_implementation'].append(background)
d['methods_ambient'][1]=d['methods_ambient'][1].replace('原检测规则对应>0；另输出估计计数≥1的资格敏感性','LIANA保留>0的原检测规则；自定义表达评分另输出估计计数≥1的资格敏感性')
p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
p=ROOT/'repro_results_static.md'
t=p.read_text(encoding='utf-8')
insert='条件置换进一步区分了这一联合证据的含义：保持另一队列观测标签时，两种选择的动态门槛同向率及相关均得到发现队列P=0.05、验证队列P=0.10。联合P=0.005由20×10的配置空间产生，而验证队列自身只有10种配置。因此，联合结果提供整体随机参照，尚不足以独立证明每个队列均存在可辨认的疾病关联。\n\n'
anchor='### 相同候选分母下的表达与原生优先级对比'
if insert not in t:t=t.replace(anchor,insert+anchor)
p.write_text(t,encoding='utf-8')
print('Integrated conditional-null and rank-universe methods/results')
p=ROOT/'repro_manuscript_en_sections.json'
d=json.loads(p.read_text(encoding='utf-8'))
d['methods_ambient']=[t.replace('LIANA retains its nonzero-count expression filter.',
    'LIANA retains the float32 normalized representation used in the raw pipeline. Positive-entry counts before and after conversion are recorded, and explicit numerical zeros are removed before its getnnz-based sparse expression proportions are calculated, enforcing a >0 detection rule.') for t in d['methods_ambient']]
p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
