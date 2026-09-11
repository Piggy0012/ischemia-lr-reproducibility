"""Retain selected previously verified figures under the new manuscript order."""
from pathlib import Path
import shutil,hashlib,json
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v2'
OLD=ROOT.parent/'outputs/revised'
mapping={'Figure6_independent_identity_review':'FigureS1_reference_identity',
         'FigureS3_doublet_score_distributions':'FigureS2_doublet_distributions',
         'FigureS4_barrier_selection_sensitivity':'FigureS3_barrier_selection',
         'Figure5_spatial_expression_context':'FigureS5_spatial_expression_context'}
audit=[]
for old,new in mapping.items():
    for ext in ['png','pdf','svg']:
        src=OLD/'figures'/f'{old}.{ext}';dst=OUT/'figures'/f'{new}.{ext}'
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
        a=hashlib.sha256(src.read_bytes()).hexdigest();b=hashlib.sha256(dst.read_bytes()).hexdigest();assert a==b
        audit.append({'source':src.relative_to(ROOT.parent).as_posix(),'output':dst.relative_to(ROOT.parent).as_posix(),'sha256':a,'unchanged':True})
(OUT/'supplement_figure_provenance.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print('Copied',len(audit),'unchanged supplementary figure files')
