"""Restore exact archived native CellChat Windows binaries without changing old libraries."""
from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'repro_v3_cellchat'
manifest=json.loads((OUT/'download_manifest.json').read_text())
for row in manifest['packages']:
    target=OUT/'downloads'/(row['package']+'_'+row['version']+'.zip')
    target.parent.mkdir(parents=True,exist_ok=True)
    if not target.exists():
        with urllib.request.urlopen(row['url'],timeout=120) as response: raw=response.read()
        assert hashlib.sha256(raw).hexdigest()==row['sha256'],f"Upstream bytes changed for {row['package']}; use the archived binary."
        target.write_bytes(raw)
    assert hashlib.sha256(target.read_bytes()).hexdigest()==row['sha256']
    print('VERIFIED',row['package'],row['version'],flush=True)
