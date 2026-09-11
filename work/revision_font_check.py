from pathlib import Path
import json
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/revised'
reports=[]
for path in (OUT/'figures').glob('*.pdf'):
    if not path.name.startswith(('Figure6_','Figure7_','FigureS3_','FigureS4_')):continue
    fonts=[]
    for page in PdfReader(path).pages:
        resources=page['/Resources'].get_object()
        for val in resources.get('/Font',{}).values():
            f=val.get_object();descendants=f.get('/DescendantFonts',[f])
            for v in descendants:
                descendant=v.get_object();descriptor=descendant.get('/FontDescriptor')
                embedded=descriptor is not None and any(k in descriptor.get_object() for k in ['/FontFile','/FontFile2','/FontFile3'])
                fonts.append({'name':str(f.get('/BaseFont')),'subtype':str(descendant.get('/Subtype')),'embedded':bool(embedded)})
    assert fonts and all(f['embedded'] for f in fonts),path.name
    reports.append({'file':path.name,'fonts':fonts})
(OUT/'figure_font_qa.json').write_text(json.dumps({'status':'PASS','checked':len(reports),'reports':reports},indent=2),encoding='utf-8')
print('Font embedding PASS',len(reports))
