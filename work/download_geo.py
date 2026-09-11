from pathlib import Path
import urllib.request, time, hashlib, json, sys, gzip

ROOT = Path(__file__).resolve().parent
def download(url, dest):
    dest = Path(dest); dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print('EXISTS', dest.name, dest.stat().st_size, flush=True); return dest
    temp = dest.with_suffix(dest.suffix+'.partial')
    print('DOWNLOAD',url,flush=True)
    req=urllib.request.Request(url,headers={'User-Agent':'Academic-public-data-audit/1.0'})
    with urllib.request.urlopen(req, timeout=180) as r, temp.open('wb') as f:
        n=0; then=time.time()
        while chunk:=r.read(1024*1024):
            f.write(chunk);n+=len(chunk)
            if time.time()-then>15: print('PROGRESS',dest.name,n,flush=True);then=time.time()
    temp.replace(dest)
    print('DONE',dest.name,dest.stat().st_size,flush=True)
    return dest

if __name__=='__main__':
    acc=sys.argv[1] if len(sys.argv)>1 else 'GSE174574'
    stem=acc[:-3]+'nnn'
    target=ROOT/'data'/acc
    urls=[f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/soft/{acc}_family.soft.gz']
    if '--matrix' in sys.argv: urls.append(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{acc}/suppl/{acc}_RAW.tar')
    records=[]
    for url in urls:
        p=download(url,target/url.rsplit('/',1)[1]);h=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
        records.append(dict(url=url,path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=h,retrieved=time.strftime('%Y-%m-%dT%H:%M:%S%z')))
        if p.name.endswith('.soft.gz'):
            txt=gzip.decompress(p.read_bytes()).decode('utf-8','replace')
            (target/(acc+'_family.soft')).write_text(txt,encoding='utf-8')
            for line in txt.splitlines():
                if line.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Sample_title','!Sample_characteristics_ch1','!Sample_source_name_ch1','!Sample_supplementary_file')): print(line,flush=True)
    (target/'download_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
