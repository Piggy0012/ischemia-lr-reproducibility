from pathlib import Path
import json,concurrent.futures
from fetch_literature import get_ref
ROOT=Path(__file__).resolve().parent
pairs=[('Dominguez2022','10.1126/science.abl5197','CellTypist'),('Yao2021','10.1016/j.cell.2021.04.021','mouse_cortex_reference'),('Yao2023','10.1038/s41586-023-06812-z','mouse_whole_brain_reference'),('Wolock2019','10.1016/j.cels.2018.11.005','Scrublet'),('Dimitrov2024','10.1038/s41556-024-01469-w','LIANA_plus'),('Efremova2020','10.1038/s41596-020-0292-x','CellPhoneDB'),('Raredon2022','10.1038/s41598-022-07959-x','Connectome'),('Hou2020','10.1038/s41467-020-18873-z','NATMI'),('Cabello2020','10.1093/nar/gkaa183','SingleCellSignalR')]
specs=[(key,doi,'https://doi.org/'+doi,role) for key,doi,role in pairs]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:refs=list(pool.map(get_ref,specs))
for r in refs:r['verified_date']='2026-09-10'
assert all(r['metadata_status']=='crossref_verified' for r in refs),refs
old=json.loads((ROOT/'references.json').read_text(encoding='utf-8'));old=[r for r in old if r['id'] not in {s[0] for s in specs}]
(ROOT/'references.json').write_text(json.dumps(old+refs,ensure_ascii=False,indent=2),encoding='utf-8')
print([(r['id'],r['title'],r['metadata_status']) for r in refs])
