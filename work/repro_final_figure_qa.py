"""Check exported dimensions and embedded fonts; visual judgments are separate."""
from pathlib import Path
import json,hashlib,sys
from PIL import Image
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/reproducibility_v2'
FIG=OUT/'figures';DATA=OUT/'figure_source_data';QA=ROOT/'repro_figure_qa'
QA.mkdir(exist_ok=True)
sys.path.insert(0,str(Path.home()/'.codex/skills/scipilot-figure-skill/scripts'))
from profile_data import profile_data
profile={'statistical_units':'11 animal libraries in paired markers; six reused representations, not independent replications',
 'chosen_chart':'paired animal points, matched-candidate correlation points and eligibility counts',
 'inference':'descriptive; no new P values or intervals for correction',
 'markers':profile_data(str(DATA/'figure4_paired_markers.tsv'),group_cols=['cell_type','gene']),
 'ranks':profile_data(str(DATA/'figure4_matched_rank_comparisons.tsv'),group_cols=['config','comparison']),
 'gates':profile_data(str(DATA/'figure4_gate_coverage.tsv'),group_cols=['config'])}
(DATA/'figure4_data_profile.json').write_text(json.dumps(profile,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
with Image.open(FIG/'Figure4_ambient_sensitivity.png') as im:im.convert('L').save(QA/'Figure4_ambient_sensitivity_grayscale.png')
reports=[]
for path in sorted(FIG.glob('*.pdf')):
    fonts=[]
    for page in PdfReader(path).pages:
        resources=page['/Resources'].get_object()
        for val in resources.get('/Font',{}).values():
            f=val.get_object()
            for v in f.get('/DescendantFonts',[f]):
                desc=v.get_object().get('/FontDescriptor')
                embedded=desc is not None and any(k in desc.get_object() for k in ['/FontFile','/FontFile2','/FontFile3'])
                fonts.append({'name':str(f.get('/BaseFont')),'embedded':bool(embedded)})
    assert fonts and all(x['embedded'] for x in fonts),path.name
    png=path.with_suffix('.png');svg=path.with_suffix('.svg');assert png.exists() and svg.exists()
    with Image.open(png) as im:size=im.size;dpi=im.info.get('dpi')
    reports.append({'figure':path.stem,'pixel_dimensions':size,'dpi':dpi,'fonts':fonts,
      'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [png,path,svg]}})
assert len(reports)==10,len(reports)
(OUT/'figure_font_and_dimension_audit.json').write_text(json.dumps({'status':'PASS','reports':reports},indent=2),encoding='utf-8')
print('All 10 PNG/PDF/SVG sets and PDF font embedding passed')
