"""Restore the full official source archive for the native CellChat binary commit."""
from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parent
meta=json.loads((ROOT/'repro_v3_cellchat_probe/corresponding_source_asset.json').read_text())
dest=ROOT.parent/meta['local_path'];dest.parent.mkdir(parents=True,exist_ok=True)
if not dest.exists():
    with urllib.request.urlopen(meta['url'],timeout=120) as response:
        payload=response.read()
    assert hashlib.sha256(payload).hexdigest()==meta['sha256'],'Archive bytes changed; use the archived release source asset'
    dest.write_bytes(payload)
assert hashlib.sha256(dest.read_bytes()).hexdigest()==meta['sha256']
print('VERIFIED official full CellChat source',meta['remote_sha'],meta['sha256'])
