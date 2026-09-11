from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import sys,hashlib,json,time
from download_geo import download
root=Path(__file__).resolve().parent
text=(root/'metadata_audit/GSE245386_family.soft').read_text(encoding='utf-8')
jobs=[];sample=None
for line in text.splitlines():
    if line.startswith('^SAMPLE = '):sample=line.split(' = ')[1]
    if line.startswith('!Sample_supplementary_file') and sample in [f'GSM784172{i}' for i in range(5)]:
        url=line.split(' = ')[1].replace('ftp://','https://');jobs.append((sample,url))
def one(job):
    sample,url=job;p=download(url,root/'data/GSE245386/raw'/url.rsplit('/',1)[1])
    return dict(sample=sample,url=url,bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest())
with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(one,jobs))
(root/'data/GSE245386/download_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
