"""Assemble v2 narrative only after required analytical result text exists."""
from pathlib import Path
import json,re
import pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2';T=OUT/'tables'

def load(name):return json.loads((ROOT/name).read_text(encoding='utf-8'))
def paragraph_list(rows):return '\n\n'.join(rows)

def single_row(frame, **selectors):
    """Fail on absent or ambiguous analytical rows instead of taking iloc[0]."""
    rows=frame
    for column,value in selectors.items():rows=rows.loc[rows[column].eq(value)]
    assert len(rows)==1,('Expected one source row',selectors,len(rows))
    return rows.iloc[0]

def ambient_methods(part, ambient, lang):
    """Use the completed model account once, followed by downstream details."""
    sections=part['methods_ambient']
    assert len(sections)>=2
    # The static first paragraph is the superseded initial-fit description.
    # Keep all subsequent normalization, detection and selection paragraphs.
    assert 'decontX 1.10.0' in sections[0] and 'maxIter=500' in sections[0]
    model=ambient['methods_'+lang].strip()
    assert all(token in model for token in ['maxIter','2000','float64'])
    if '[@Yang2020]' not in model:model+=' [@Yang2020]'
    return model+'\n\n'+paragraph_list(sections[1:])

def insert_before_once(text, anchor, addition):
    """Keep completed analytical narratives visible through later insertions."""
    assert addition.strip() and addition not in text
    assert text.count(anchor)==1,('Expected one narrative anchor',anchor)
    text=text.replace(anchor,addition+'\n\n'+anchor,1)
    assert text.count(addition)==1
    return text

def tables(lang):
    if lang=='zh':
        table1='表1　数据来源与在本研究中的用途\n\n| GEO | 对比及样本单位 | Sham / MCAO | 质控数 | 用途 |\n| --- | --- | --- | --- | --- |\n'
        rows=[['GSE174574','急性缺血，动物文库','3 / 3','58,333细胞','主要队列'],['GSE245386','野生型，动物文库','3 / 2','57,022细胞','原外部队列'],['GSE332910','同侧纹状体单核，作者报告生物学重复','3 / 3','116,308核','新增描述性扩展'],['GSE163752','缺血侧与对侧，独立样本池','6对池','不适用','分选效应与区间'],['GSE233814','Control及D1/D3/D7切片记录','5张切片','11,969空间点','单张补充表达图']]
        head='表2　相同目标候选下的表达与优先级比较\n\n| 选择 | 指标 | N | 全部同向 | 非零同向 | ρ |\n| --- | --- | --- | --- | --- | --- |\n'
        names={'primary':'原规则','reference_singlet':'参考支持'}
        metric_names={'custom_score':'表达共可用性','native_magnitude_priority':'原生LIANA 1.10.0','diagnostic_unique_column_priority':'独特分数列诊断','fixed':'固定全局母集诊断'}
        foot='N为固定目标候选数。全部同向以N为分母；非零同向仅以两队列效应均非零者为分母。ρ为包括全部候选的Spearman相关。这些是同一批动物及候选的表示比较，无新增方法间显著性结论。'
    else:
        table1='Table 1. Data sources and analytical roles\n\n| GEO | Contrast and sample unit | Sham / MCAO | QC retained | Role |\n| --- | --- | --- | --- | --- |\n'
        rows=[['GSE174574','Acute ischemia; animal libraries','3 / 3','58,333 cells','Primary cohort'],['GSE245386','Wild type; animal libraries','3 / 2','57,022 cells','Original validation'],['GSE332910','Ipsilateral striatum; author-reported biological replicates','3 / 3','116,308 nuclei','Descriptive extension'],['GSE163752','Ipsilateral vs contralateral; independent pools','6 paired pools','Not applicable','Sorted effect intervals'],['GSE233814','Control and D1/D3/D7 section records','5 sections','11,969 spots','One supplementary map']]
        head='Table 2. Expression and priority agreement on matched targets\n\n| Selection | Metric | N | All concordant | Nonzero concordant | ρ |\n| --- | --- | --- | --- | --- | --- |\n'
        names={'primary':'Primary','reference_singlet':'Reference'}
        metric_names={'custom_score':'Expression coavailability','native_magnitude_priority':'Native LIANA 1.10.0','diagnostic_unique_column_priority':'Unique-column diagnostic','fixed':'Fixed-global diagnostic'}
        foot='N denotes the fixed target-candidate count. All-candidate concordance uses N; nonzero concordance uses effects nonzero in both cohorts. Spearman correlation includes all N candidates. These representations reuse the same animals and targets and do not provide new formal between-method significance evidence.'
    table1+='\n'.join('| '+' | '.join(row)+' |' for row in rows)
    null=pd.read_csv(T/'null_fixed_common_summary.tsv',sep='\t',float_precision='round_trip')
    allocations=pd.read_csv(T/'null_fixed_common_all_allocations.tsv.gz',sep='\t',float_precision='round_trip')
    diag=pd.read_csv(T/'rank_diagnostic_raw_summary.tsv',sep='\t',float_precision='round_trip')
    bg=pd.read_csv(T/'rank_background_summary.tsv',sep='\t',float_precision='round_trip')
    for config in names:
        common_n=None
        for metric in metric_names:
            if metric=='custom_score':
                row=single_row(allocations,config=config,metric=metric,is_observed=True)
            elif metric=='fixed':
                row=single_row(bg,config=config,metric='unique_fixed_global_priority')
            else:
                row=single_row(diag,count_source='raw',config=config,metric=metric)
            counts=[row[column] for column in ['n_candidates','n_same_direction','n_nonzero_both']]
            assert all(pd.notna(value) and float(value).is_integer() for value in counts)
            n,same,nz=map(int,counts);rho=float(row.spearman_rho)
            assert 0<=same<=nz<=n and nz>0 and -1<=rho<=1,(config,metric,n,same,nz,rho)
            assert abs(float(row.same_direction_all_fraction)-same/n)<1e-12
            assert abs(float(row.same_direction_nonzero_fraction)-same/nz)<1e-12
            if common_n is None:common_n=n
            assert n==common_n,('Unmatched Table 2 candidate counts',config,metric,n,common_n)
            # Also reconcile custom/native rows against the inferential export.
            if metric in ['custom_score','native_magnitude_priority']:
                null_metric='magnitude_priority' if metric=='native_magnitude_priority' else metric
                for statistic,value in [('spearman_rho',rho),('same_direction_all_fraction',same/n),
                                        ('same_direction_nonzero_fraction',same/nz)]:
                    check=single_row(null,config=config,metric=null_metric,statistic=statistic,
                                     test_family='A_fixed_common_metrics')
                    assert check.observed_n_candidates==n
                    assert abs(float(check.observed)-value)<1e-12,(config,metric,statistic)
            head+='| '+ ' | '.join([names[config],metric_names[metric],str(n),f'{same}/{n} ({same/n:.1%})',f'{same}/{nz} ({same/nz:.1%})',f'{rho:.3f}'])+' |\n'
    return table1+'\n\n'+head+'\n'+foot

def main():
    front=load('repro_manuscript_frontmatter.json');caps=load('repro_caption_static.json');author=load('repro_author.json')
    ambient=load('repro_ambient_narrative.json')
    third=load('repro_third_liana_narrative.json')
    # These files must be produced from complete actual runs, never promises.
    assert (T/'ambient_rank_matched_audit.json').exists()
    assert ambient['status']=='complete' and third['status']=='complete'
    for lang in ['zh','en']:
        part=load('repro_main_sections.json' if lang=='zh' else 'repro_manuscript_en_sections.json')
        name=author['name'] if lang=='zh' else author['preferred_romanization']
        institution=author['affiliation_zh']+'（'+author['affiliation']+'）' if lang=='zh' else author['affiliation']
        text=f'# {front["title_"+lang]}\n\n{name}\n\n{institution}\n\n'
        text+=('通讯邮箱：' if lang=='zh' else 'Correspondence: ')+author['corresponding_email']+'\n\n'
        text+='[ORCID 0009-0005-9705-0585](https://orcid.org/0009-0005-9705-0585)\n\n'
        text+=('## 摘要' if lang=='zh' else '## Abstract')+'\n\n'+front['abstract_'+lang]+'\n\n'
        text+=('关键词：' if lang=='zh' else 'Keywords: ')+front['keywords_'+lang]+'\n\n'
        text+=('## 引言' if lang=='zh' else '## Introduction')+'\n\n'+paragraph_list(part['introduction'])+'\n\n'
        text+=('## 方法' if lang=='zh' else '## Methods')+'\n\n'
        headers_zh=['数据和细胞选择','动物级评分和共同候选','保留依赖结构的标签零基线','聚合实现和排名母集诊断','环境RNA与数值收敛敏感性','分选和空间补充分析']
        headers_en=['Data and cell selection','Animal-level scores and matched candidates','Dependence-preserving label nulls','Aggregation implementation and ranking universes','Ambient RNA and numerical convergence','Sorted and spatial supplementary analyses']
        keys=['methods_design','methods_scores','methods_null','methods_implementation','methods_ambient','methods_auxiliary']
        for key,h in zip(keys,headers_zh if lang=='zh' else headers_en):
            body=ambient_methods(part,ambient,lang) if key=='methods_ambient' else paragraph_list(part[key])
            text+='### '+h+'\n\n'+body+'\n\n'
            if key=='methods_scores':text+=third['methods_'+lang]+'\n\n'
        assert text.count(third['methods_'+lang])==1
        if lang=='zh':
            result=(ROOT/'repro_results_static.md').read_text(encoding='utf-8')
            result=insert_before_once(result,'### 第三队列显示',ambient['results_zh'])
            result=insert_before_once(result,'### 分选效应',third['results_zh'])
            discussion=(ROOT/'repro_discussion_zh.md').read_text(encoding='utf-8')
            discussion=insert_before_once(discussion,'本研究的限制集中',ambient['discussion_zh'])
        else:
            core=(ROOT/'repro_manuscript_en_core.md').read_text(encoding='utf-8')
            result=core.split('## Results\n',1)[1].split('\n## Discussion',1)[0].strip()
            discussion=core.split('\n## Discussion\n',1)[1].strip()
            result=insert_before_once(result,'### A third cohort',ambient['results_en'])
            result=insert_before_once(result,'Sorted-data summaries',third['results_en'])
            discussion=insert_before_once(discussion,'Several limitations delimit',ambient['discussion_en'])
        assert result.count(third['results_'+lang])==1 and result.count(ambient['results_'+lang])==1
        text+=('## 结果' if lang=='zh' else '## Results')+'\n\n'+result+'\n\n'+tables(lang)+'\n\n'
        text+=('## 讨论' if lang=='zh' else '## Discussion')+'\n\n'+discussion+'\n\n'
        text+=('## 结论' if lang=='zh' else '## Conclusion')+'\n\n'+front['conclusion_'+lang]+'\n\n'
        text+=('## 数据获取' if lang=='zh' else '## Data availability')+'\n\n'
        links=', '.join(f'[{acc}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc})' for acc in ['GSE174574','GSE245386','GSE332910','GSE163752','GSE233814'])
        text+=('原始公开数据可从GEO获取：' if lang=='zh' else 'Public source data are available from GEO: ')+links+'.\n\n'
        record=OUT/'archive_record.json'
        if record.exists():
            d=json.loads(record.read_text(encoding='utf-8'));assert d['public_verified'] is True
            assert re.fullmatch(r'10\.5281/zenodo\.\d+',d['doi'])
            text+=('## 代码与派生数据\n\n代码、环境版本、细胞汇总计数及全部检验结果已归档于' if lang=='zh' else '## Code and derived data\n\nCode, environment versions, cell-summary counts, and complete test results are archived at ')+f'[Zenodo](https://doi.org/{d["doi"]}).\n\n'
        text+=('## 伦理说明\n\n本研究仅重分析公开动物数据，未开展新的动物或人体实验。原研究的实验审批与许可见相应来源文献。' if lang=='zh' else '## Ethics statement\n\nThis study reanalyzed public animal data and performed no new animal or human experiments. Approval and consent for source experiments are reported by the original studies.')+'\n\n'
        text+=('## 图注' if lang=='zh' else '## Figure legends')+'\n\n'
        legends={**caps[lang],'Figure4':ambient['caption_'+lang],'Figure5':third['caption_'+lang]}
        for label in ['Figure1','Figure2','Figure3','Figure4','Figure5','FigureS1','FigureS2','FigureS3','FigureS4','FigureS5']:
            text+=legends[label]+'\n\n'
        text+=('## 参考文献' if lang=='zh' else '## References')+'\n\n<!-- REFERENCES_GENERATED -->\n'
        (ROOT/f'repro_manuscript_{lang}_source.md').write_text(text,encoding='utf-8')
    (OUT/'manuscript_readiness.json').write_text(json.dumps({'analysis_text_complete':True,
        'public_archive_verified':(OUT/'archive_record.json').exists(),
        'author_affiliation_and_correspondence_supplied':True,
        'submission_ready':False,'reason':'Public deposition and author/journal submission metadata must be verified separately.'},indent=2),encoding='utf-8')
    print('Prepared Chinese and English v2 manuscript sources')

if __name__=='__main__':main()
