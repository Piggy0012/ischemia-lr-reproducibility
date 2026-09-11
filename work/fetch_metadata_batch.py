from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request,gzip,time
root=Path(__file__).resolve().parent/'metadata_audit';root.mkdir(exist_ok=True)
accs=['GSE233812','GSE233813','GSE233815','GSE154396','GSE245386','GSE234052']
def one(acc):
    url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz'
    try:
        data=urllib.request.urlopen(url,timeout=90).read()
        (root/(acc+'_family.soft.gz')).write_bytes(data)
        txt=gzip.decompress(data).decode('utf-8','replace');(root/(acc+'_family.soft')).write_text(txt,encoding='utf-8')
        print(acc,'OK',len(data),flush=True)
    except Exception as e: print(acc,'ERROR',str(e),flush=True)
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(one,accs))
try:
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=34496660,35688114,37183260,37063847,39953063,34836962&retmode=xml'
    (root/'pubmed_metadata.xml').write_bytes(urllib.request.urlopen(url,timeout=90).read());print('PUBMED OK',flush=True)
except Exception as e:print('PUBMED ERROR',e,flush=True)
