"""Fetch canonical bibliographic metadata for the methodological revision."""
from pathlib import Path
import json,urllib.parse,urllib.request,time
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'repro_literature';OUT.mkdir(exist_ok=True)
NEW={'Yang2020':'10.1186/s13059-020-1950-6','Kolde2012':'10.1093/bioinformatics/btr709',
     'Cesaro2025':'10.1093/nargab/lqaf084','Ku2026':'10.1186/s13059-026-04063-5',
     'Zhang2026':'10.1002/advs.77547','Hou2025':'10.1038/s41467-025-66272-z'}
refs=json.loads((ROOT/'literature/references.json').read_text(encoding='utf-8'))
for key,doi in NEW.items():
    path=OUT/(key+'_crossref.json')
    url='https://api.crossref.org/works/'+urllib.parse.quote(doi,safe='')
    if not path.exists():
        req=urllib.request.Request(url,headers={'User-Agent':'Research-reproducibility-audit/2.0'})
        with urllib.request.urlopen(req,timeout=60) as r:raw=r.read()
        d=json.loads(raw);assert d['status']=='ok'
        path.write_bytes(raw)
    d=json.loads(path.read_bytes())['message'];assert d['DOI'].lower()==doi.lower()
    date=d.get('published',d.get('published-online',{})).get('date-parts',[[None]])[0]
    rec={'id':key,'title':d['title'][0],
         'authors':[' '.join(filter(None,[a.get('given'),a.get('family')])) for a in d.get('author',[])],
         'year':date[0],'journal':d.get('container-title',[''])[0],
         'volume':d.get('volume'),'issue':d.get('issue'),'pages':d.get('page'),
         'article_number':d.get('article-number'),'doi':doi,'url':'https://doi.org/'+doi,
         'metadata_source':url,'metadata_status':'crossref_verified','verified_date':'2026-09-11'}
    refs=[r for r in refs if r['id']!=key]+[rec]
    print(key,rec['year'],rec['title'])
(OUT/'references.json').write_text(json.dumps(refs,indent=2,ensure_ascii=False),encoding='utf-8')
