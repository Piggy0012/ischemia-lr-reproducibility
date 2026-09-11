"""Record the root review actually performed on rendered manuscript pages."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/reproducibility_v2'
final=ROOT/'repro_docx_release_zh';reviewed=ROOT/'repro_docx_final_zh';old=ROOT/'repro_docx_render_zh'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=[]
for i in range(1,25):
    p=final/f'page-{i}.png';prior=old/p.name
    assert sha(p)==sha(reviewed/p.name),'Release layout changed after visual review'
    identical=sha(p)==sha(prior)
    # All old pages 1–20 were displayed individually in the root review.
    # Final changed pages 2/13/14 and final pages 21–24 were displayed again.
    direct=i in [2,13,14,21,22,23,24]
    assert direct or (i<=20 and identical)
    rows.append({'page':i,'png_sha256':sha(p),'status':'PASS',
      'basis':'visually inspected final-layout page; release PNG byte-identical' if direct else 'visually inspected prior-render page; release PNG byte-identical',
      'checks':['legible text and symbols','no clipping or overlap','table/header or figure placement','page number']})
report={'status':'PASS','language':'zh','n_pages':24,'reviewed_utc':datetime.now(timezone.utc).isoformat(),
 'docx_sha256':sha(OUT/'manuscript_zh.docx'),'rendered_pdf_sha256':sha(final/'manuscript_zh.pdf'),
 'pages':rows,'observations':['Table 1 spans pages 8–9 with a repeated header and intact rows.',
 'Ten figures are presented individually on pages 15–24.',
 'Reference entries remain intact across pages 13–14.'],
 'scope':'rendered layout only; no claim of submission or public deposition'}
(OUT/'validation/document_visual_review_zh.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('Recorded actual 24-page Chinese visual review')
