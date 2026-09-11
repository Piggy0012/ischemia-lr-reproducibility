import requests,json,hashlib
from pathlib import Path
p=Path('work/identity_reference_review');m=requests.get('https://celltypist.cog.sanger.ac.uk/models/models.json',timeout=60).json();(p/'models.json').write_text(json.dumps(m,indent=2))
rows=[]
for r in m['models']:
 if r['filename'] not in ['Mouse_Isocortex_Hippocampus.pkl','Mouse_Whole_Brain.pkl']:continue
 dst=p/r['filename']
 if not dst.exists():
  res=requests.get(r['url'],timeout=90);res.raise_for_status();dst.write_bytes(res.content)
 r.update(bytes=dst.stat().st_size,sha256=hashlib.sha256(dst.read_bytes()).hexdigest());rows.append(r);print(json.dumps(r),flush=True)
(p/'manifest.json').write_text(json.dumps(rows,indent=2))
