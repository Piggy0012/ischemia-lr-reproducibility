from pathlib import Path
import concurrent.futures, json, urllib.request, urllib.parse, datetime, hashlib

ROOT=Path(__file__).resolve().parent
SPECS=[
('Liu2026','10.1186/s13062-026-00912-2','https://link.springer.com/article/10.1186/s13062-026-00912-2','closest_comparator'),
('Hill2025','10.1038/s41467-025-65487-4','https://www.nature.com/articles/s41467-025-65487-4','endfoot_localization'),
('Spitzer2022','10.1007/s00401-022-02452-1','https://pubmed.ncbi.nlm.nih.gov/35752654/','GSE163752_source_SPP1_prior'),
('Zheng2022','10.1177/0271678X211026770','https://pubmed.ncbi.nlm.nih.gov/34496660/','GSE174574_source'),
('Munji2019','10.1038/s41593-019-0497-x','https://www.nature.com/articles/s41593-019-0497-x','BBB_signatures'),
('Vanlandewijck2018','10.1038/nature25739','https://www.nature.com/articles/nature25739','vascular_annotation'),
('Armulik2010','10.1038/nature09522','https://www.nature.com/articles/nature09522','pericyte_BBB_function'),
('Nitta2003','10.1083/jcb.200302070','https://pubmed.ncbi.nlm.nih.gov/12743111/','Cldn5_function'),
('BenZvi2014','10.1038/nature13324','https://pmc.ncbi.nlm.nih.gov/articles/PMC4134871/','Mfsd2a_function'),
('Andreone2017','10.1016/j.neuron.2017.03.043','https://pubmed.ncbi.nlm.nih.gov/28416077/','transcytosis_function'),
('DelToro2010','10.1182/blood-2010-02-270819','https://pubmed.ncbi.nlm.nih.gov/20705756/','tip_cell_genes'),
('Jin2021','10.1038/s41467-021-21246-9','https://www.nature.com/articles/s41467-021-21246-9','CellChat_method'),
('Dimitrov2022','10.1038/s41467-022-30755-0','https://www.nature.com/articles/s41467-022-30755-0','LR_benchmark_LIANA'),
('Browaeys2020','10.1038/s41592-019-0667-5','https://www.nature.com/articles/s41592-019-0667-5','NicheNet_method'),
('Squair2021','10.1038/s41467-021-25960-2','https://www.nature.com/articles/s41467-021-25960-2','pseudobulk_false_discoveries'),
('Crowell2020','10.1038/s41467-020-19894-4','https://pubmed.ncbi.nlm.nih.gov/33257685/','muscat_multi_sample'),
('Shao2021','10.1093/bib/bbaa269','https://academic.oup.com/bib/article/22/4/bbaa269/5955941','CellTalkDB_resource'),
('Quenault2017','10.1093/brain/aww260','https://academic.oup.com/brain/article/140/1/146/2528261','endothelial_activation_genes'),
('Xu2026','10.1002/cns.70886','https://pmc.ncbi.nlm.nih.gov/articles/PMC13122273/','SPP1_prior_GVU'),
]

def fetch(url):
    req=urllib.request.Request(url,headers={'User-Agent':'Codex-literature-audit/1.0 (public research reproducibility)'})
    with urllib.request.urlopen(req,timeout=45) as r: return r.read()

def get_ref(spec):
    key,doi,source,role=spec
    api='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
    try:
        raw=fetch(api); obj=json.loads(raw)['message']
        (ROOT/(key+'_crossref.json')).write_bytes(raw)
        date=obj.get('published-print',obj.get('published',{})).get('date-parts',[[]])[0]
        return {'id':key,'title':obj.get('title',[''])[0],
            'authors':[(' '.join([a.get('given',''),a.get('family','')])).strip() for a in obj.get('author',[])],
            'year':date[0] if date else None,'journal':obj.get('container-title',[''])[0],
            'volume':obj.get('volume'),'issue':obj.get('issue'),'pages':obj.get('page'),'article_number':obj.get('article-number'),
            'doi':doi,'url':source,'role':role,'verified_source':source,
            'metadata_source':api,'verified_date':'2026-09-09','metadata_status':'crossref_verified'}
    except Exception as e:
        return {'id':key,'doi':doi,'url':source,'role':role,'verified_source':source,'metadata_status':'fetch_failed','error':str(e)}

if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        refs=list(pool.map(get_ref,SPECS))
    (ROOT/'references.json').write_text(json.dumps(refs,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'references':len(refs),'failed':[r['id'] for r in refs if r['metadata_status']!='crossref_verified']}))
    for pmc in ['PMC6858546','PMC4314527']:
        url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML'
        try:
            data=fetch(url); (ROOT/(pmc+'.xml')).write_bytes(data); print(pmc,len(data))
        except Exception as e: print(pmc,str(e))
    try:
        commit=json.loads(fetch('https://api.github.com/repos/saezlab/liana-py/commits?path=src/liana/resource/omni_resource.csv&per_page=1'))[0]
        version={'upstream_commit':commit['sha'],'commit_url':commit['html_url']}
    except Exception as e: version={'commit_lookup_error':str(e)}
    version.update({'retrieved':'2026-09-09','upstream_url':'https://raw.githubusercontent.com/saezlab/liana-py/main/src/liana/resource/omni_resource.csv','mouse_rows':3989,
        'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'omni_resource.csv',ROOT/'mouseconsensus.csv']}})
    (ROOT/'lr_resource_manifest.json').write_text(json.dumps(version,indent=2),encoding='utf-8')
