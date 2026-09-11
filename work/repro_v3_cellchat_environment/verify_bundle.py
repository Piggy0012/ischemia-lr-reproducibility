"""Verify an extracted environment archive using only the Python standard library."""
from pathlib import Path
import argparse,csv,hashlib,json,zipfile
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main(root):
    root=root.resolve(); env=root/'environment'
    m=json.loads((env/'environment_manifest.json').read_text(encoding='utf-8'))
    for r in m['files']:
        p=(env/r['path']).resolve();assert p.is_relative_to(env)
        assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],r['path']
    with (env/'offline_binary_manifest.tsv').open(encoding='utf-8') as f:bins=list(csv.DictReader(f,delimiter='\t'))
    assert len(bins)==42
    for r in bins:
        p=(root/r['archive_path']).resolve();assert p.is_relative_to(root)
        assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['Package']
        with zipfile.ZipFile(p) as z:assert z.testzip() is None,r['Package']
    print(json.dumps({'status':'passed','metadata_files':len(m['files']),'binary_packages':len(bins),'installation_performed':False}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('extracted_bundle',type=Path);a=p.parse_args();main(a.extracted_bundle)
