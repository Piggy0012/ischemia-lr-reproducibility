from pathlib import Path
import gzip,json,hashlib
from concurrent.futures import ThreadPoolExecutor
from download_geo import download
root=Path(__file__).resolve().parent
soft=(root/'metadata_audit/GSE233815_family.soft').read_text(encoding='utf-8')
jobs=[];sample=None
for line in soft.splitlines():
    if line.startswith('^SAMPLE = '):sample=line.split(' = ')[1]
    if line.startswith('!Sample_supplementary_file') and sample in ['GSM7437221','GSM7437222','GSM7437223','GSM7437224','GSM7437225']:
        url=line.split(' = ')[1].replace('ftp://','https://')
        if not any(x in url.lower() for x in ['.tif','.jpg','.png']):jobs.append((sample,url))
def one(job):
    sample,url=job;p=download(url,root/'data/GSE233814/raw'/url.rsplit('/',1)[1])
    return dict(sample=sample,url=url,bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest())
with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(one,jobs))
(root/'data/GSE233814/download_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
for p in (root/'data/GSE233814/raw').glob('*.json.gz'):
    obj=json.loads(gzip.decompress(p.read_bytes()));print(p.name,list(obj)[:15],str(obj)[:1000],flush=True)
