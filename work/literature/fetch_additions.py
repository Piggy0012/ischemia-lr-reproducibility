from pathlib import Path
import json,concurrent.futures
from fetch_literature import get_ref,fetch
ROOT=Path(__file__).resolve().parent
specs=[('Ruan2023','10.1186/s12974-023-02941-4','https://pmc.ncbi.nlm.nih.gov/articles/PMC10687904/','GSE245386_source'),('RuanCorrection2025','10.1186/s12974-025-03610-4','https://doi.org/10.1186/s12974-025-03610-4','GSE245386_correction'),('Zucha2024','10.1073/pnas.2404203121','https://doi.org/10.1073/pnas.2404203121','GSE233814_source'),('Muzellec2023','10.1093/bioinformatics/btad547','https://pmc.ncbi.nlm.nih.gov/articles/PMC10502239/','PyDESeq2'),('Arbaizar2023','10.1186/s12974-023-02888-6','https://pmc.ncbi.nlm.nih.gov/articles/PMC10494365/','endothelial_activation'),('Qi2003','10.1038/nm846','https://www.nature.com/articles/nm846','TIMP3_KDR_primary_evidence'),('Yepes2003','10.1172/JCI19212','https://pmc.ncbi.nlm.nih.gov/articles/PMC259131/','tPA_LRP_primary_evidence'),('Mikelis2009','10.1096/fj.08-117564','https://faseb.onlinelibrary.wiley.com/doi/10.1096/fj.08-117564','PTN_PTPRZ1_not_PTPRB'),('Love2014','10.1186/s13059-014-0550-8','https://doi.org/10.1186/s13059-014-0550-8','DESeq2')]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:refs=list(pool.map(get_ref,specs))
old=json.loads((ROOT/'references.json').read_text(encoding='utf-8'));old=[r for r in old if r['id'] not in {s[0] for s in specs}]
(ROOT/'references.json').write_text(json.dumps(old+refs,ensure_ascii=False,indent=2),encoding='utf-8');print([(r['id'],r['metadata_status']) for r in refs],flush=True)
for pmc in ['PMC8721774','PMC10687904','PMC9288377','PMC10494365','PMC259131']:
 p=ROOT/(pmc+'.xml')
 if not p.exists():
  try:p.write_bytes(fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML'));print(pmc,p.stat().st_size,flush=True)
  except Exception as e:print(pmc,str(e),flush=True)
