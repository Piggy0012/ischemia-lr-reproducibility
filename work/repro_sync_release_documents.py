"""Sync final reviewed documents without replacing relocated analytical tables."""
from pathlib import Path
import json,shutil,hashlib
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/reproducibility_v2'
STAGE=ROOT/'repro_delivery/analysis'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
scientific=json.loads((OUT/'validation/manuscript_scientific_review.json').read_text(encoding='utf-8'))
copies=[]
for lang in ['zh','en']:
    v=json.loads((OUT/f'validation/document_visual_review_{lang}.json').read_text(encoding='utf-8'))
    assert str(v['status']).upper()=='PASS'
    # Reviewers record the exact delivered DOCX rather than an earlier draft.
    text=json.dumps(v)
    assert sha(OUT/f'manuscript_{lang}.docx') in text
for pattern in ['manuscript_*.docx','manuscript_*.md','references_*.json','document_build_audit.json',
                'README.md','revision_response_12_points_zh.md','manuscript_readiness.json','figure_font_and_dimension_audit.json']:
    for src in OUT.glob(pattern):copies.append((src,STAGE/'outputs/reproducibility_v2'/src.name))
for name in ['document_visual_review_zh.json','document_visual_review_en.json','manuscript_scientific_review.json']:
    copies.append((OUT/'validation'/name,STAGE/'outputs/reproducibility_v2/validation'/name))
for name in ['repro_build_docx.py','repro_prepare_manuscripts.py','repro_finalize_narrative.py',
             'repro_record_zh_visual_review.py','repro_final_figure_qa.py','repro_sync_release_documents.py',
             'repro_manuscript_zh_source.md','repro_manuscript_en_source.md','repro_manuscript_frontmatter.json',
             'repro_main_sections.json','repro_manuscript_en_sections.json','repro_caption_static.json']:
    copies.append((ROOT/name,STAGE/'work'/name))
report=[]
for src,dest in copies:
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    assert sha(src)==sha(dest)
    report.append({'source':str(src.relative_to(ROOT.parent)),'staged':str(dest.relative_to(STAGE)),
                   'bytes':src.stat().st_size,'sha256':sha(src)})
payload={'status':'PASS','scope':'final reviewed documents and presentation metadata only; relocated tables retained',
         'files':report,'public_deposition_complete':False}
for folder in [OUT/'validation',STAGE/'outputs/reproducibility_v2/validation']:
    (folder/'final_document_sync.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
print('Synchronized',len(report),'reviewed document and source files')
