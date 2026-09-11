"""Build bilingual manuscript DOCX from audited v2 Markdown sources.

Use the bundled document runtime. Render and inspect both documents separately
after this authoring step; successful OOXML serialization is not visual QA.
"""
from pathlib import Path
import re,json,html
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from PIL import Image
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs/reproducibility_v3'
REFS={r['id']:r for r in json.loads((ROOT/'repro_literature/references.json').read_text(encoding='utf-8'))}
EXTRA_REFS=ROOT/'repro_v3_cellchat_references.json'
if EXTRA_REFS.exists():
    for entry in json.loads(EXTRA_REFS.read_text(encoding='utf-8')):
        assert entry['id'] not in REFS, entry['id']
        REFS[entry['id']]=entry

def clean(s):return ' '.join(html.unescape(re.sub(r'<[^>]+>','',str(s))).split()).replace('α v β 3','αvβ3')

def build(lang, kind):
    raw=(ROOT/f'repro_v3_{kind}_{lang}_source.md').read_text(encoding='utf-8')
    assert '<!-- PENDING' not in raw and '本稿不宣称' not in raw
    order=[]
    def cite(m):
        key=m.group(1);assert key in REFS and REFS[key]['metadata_status'] in {'crossref_verified','primary_journal_verified'},key
        if key not in order:order.append(key)
        return f'[{order.index(key)+1}](https://doi.org/{REFS[key]["doi"]})'
    md=re.sub(r'\[@(\w+)\]',cite,raw)
    bib=[]
    for i,key in enumerate(order,1):
        r=REFS[key];authors=', '.join(r['authors'][:6])+(', et al.' if len(r['authors'])>6 else '.')
        locator=r.get('pages') or r.get('article_number') or ''
        vi=(str(r.get('volume') or '')+(':'+str(locator) if locator else '')).strip(':')
        bib.append(f'{i}. {clean(authors)} {clean(r["title"])}. {clean(r["journal"])}. {r["year"]}'+(f';{vi}' if vi else '')+f'. [doi:{r["doi"]}](https://doi.org/{r["doi"]})')
    assert md.count('<!-- REFERENCES_GENERATED -->')==1
    md=md.replace('<!-- REFERENCES_GENERATED -->','\n\n'.join(bib))
    (OUT/f'{kind}_{lang}.md').write_text(md,encoding='utf-8')
    (OUT/f'{kind}_references_{lang}.json').write_text(json.dumps([REFS[k] for k in order],ensure_ascii=False,indent=2),encoding='utf-8')
    doc=Document()
    for e in doc.styles.element.xpath('.//w:pBdr'):e.getparent().remove(e)
    sec=doc.sections[0];sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    sec.top_margin=sec.bottom_margin=Inches(.72);sec.left_margin=sec.right_margin=Inches(.65)
    for name in ['Normal','Title','Subtitle','Heading 1','Heading 2','Heading 3','Caption']:
        st=doc.styles[name];st.font.name='Times New Roman';st.font.color.rgb=RGBColor(0,0,0)
        st._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'SimSun')
    normal=doc.styles['Normal'];normal.font.size=Pt(11)
    normal.paragraph_format.line_spacing=1.3;normal.paragraph_format.space_after=Pt(6)
    normal.paragraph_format.widow_control=True
    for name,size in [('Title',17),('Heading 1',14),('Heading 2',12),('Heading 3',11)]:
        st=doc.styles[name];st.font.size=Pt(size);st.font.bold=True
        st.paragraph_format.space_before=Pt(12);st.paragraph_format.space_after=Pt(6)
        st.paragraph_format.keep_with_next=True
        st._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'SimHei')
    doc.styles['Title'].paragraph_format.space_before=Pt(0)
    footer=sec.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
    fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE')
    run=OxmlElement('w:r');txt=OxmlElement('w:t');txt.text='1';run.append(txt);fld.append(run);footer._p.append(fld)

    def inline(p,text):
        for part in re.split(r'(\*\*.*?\*\*|\[[^\]]+\]\([^)]+\)|`[^`]+`)',text):
            if not part:continue
            if part.startswith('**') and part.endswith('**'):p.add_run(part[2:-2]).bold=True
            elif part.startswith('`') and part.endswith('`'):p.add_run(part[1:-1])
            elif re.fullmatch(r'\[[^\]]+\]\([^)]+\)',part):
                label,url=re.match(r'\[([^\]]+)\]\(([^)]+)\)',part).groups()
                h=OxmlElement('w:hyperlink');h.set(qn('r:id'),p.part.relate_to(url,RT.HYPERLINK,is_external=True))
                r=OxmlElement('w:r');pr=OxmlElement('w:rPr');c=OxmlElement('w:color');c.set(qn('w:val'),'000000');pr.append(c);r.append(pr)
                t=OxmlElement('w:t');t.text='['+label+']' if label.isdigit() else label;r.append(t);h.append(r);p._p.append(h)
            else:p.add_run(part)

    def table(lines):
        rows=[[c.strip() for c in l.strip().strip('|').split('|')] for l in lines if not re.match(r'^\|[\s:|\-]+\|$',l)]
        n=len(rows[0]);assert all(len(r)==n for r in rows)
        widths={4:[1.25,2.4,1.4,2.15],5:[1.05,2.1,1.1,1.05,1.9],6:[.85,2.05,.55,1.25,1.7,.8]}.get(n)
        assert widths and abs(sum(widths)-7.2)<.001
        tb=doc.add_table(rows=1,cols=n);tb.alignment=WD_TABLE_ALIGNMENT.CENTER;tb.autofit=False
        for i,col in enumerate(tb.columns):col.width=Inches(widths[i])
        for ri,row in enumerate(rows):
            cells=tb.rows[0].cells if ri==0 else tb.add_row().cells
            for ci,text in enumerate(row):
                cell=cells[ci];cell.width=Inches(widths[ci]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
                p=cell.paragraphs[0];p.paragraph_format.line_spacing=1.1
                p.paragraph_format.space_before=Pt(4);p.paragraph_format.space_after=Pt(4)
                if ci>1:p.alignment=WD_ALIGN_PARAGRAPH.CENTER
                inline(p,text)
                for r in p.runs:r.font.size=Pt(9);r.bold=ri==0
                tcpr=cell._tc.get_or_add_tcPr();borders=OxmlElement('w:tcBorders')
                for side in ['top','left','bottom','right']:
                    el=OxmlElement('w:'+side);el.set(qn('w:val'),'single');el.set(qn('w:sz'),'4');el.set(qn('w:color'),'D9D9D9');borders.append(el)
                tcpr.append(borders)
                margins=OxmlElement('w:tcMar')
                for side in ['top','left','bottom','right']:
                    el=OxmlElement('w:'+side);el.set(qn('w:w'),'65');el.set(qn('w:type'),'dxa');margins.append(el)
                tcpr.append(margins)
                if ri==0:
                    shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'E6EEF3');tcpr.append(shade)
            prop=tb.rows[ri]._tr.get_or_add_trPr();prop.append(OxmlElement('w:cantSplit'))
            if ri==0:prop.append(OxmlElement('w:tblHeader'))
        doc.add_paragraph()

    lines=md.splitlines();i=0;in_references=False
    while i<len(lines):
        line=lines[i].strip();i+=1
        if not line:continue
        if line.startswith('<!--'):continue
        if line.startswith('|'):
            block=[line]
            while i<len(lines) and lines[i].strip().startswith('|'):block.append(lines[i]);i+=1
            table(block);continue
        if line.startswith('# '):inline(doc.add_paragraph(style='Title'),line[2:].replace('与优先级','\n与优先级').replace('ligand–receptor','ligand\u2060–\u2060receptor'));continue
        if line.startswith('## '):
            in_references=line[3:] in ['References','参考文献']
            inline(doc.add_paragraph(style='Heading 1'),line[3:]);continue
        if line.startswith('### '):inline(doc.add_paragraph(style='Heading 2'),line[4:]);continue
        if line.startswith('- '):inline(doc.add_paragraph(style='List Bullet'),line[2:]);continue
        if re.match(r'^\d+\. ',line) and not in_references:
            inline(doc.add_paragraph(style='List Number'),re.sub(r'^\d+\. ', '', line));continue
        p=doc.add_paragraph();inline(p,line)
        if re.match(r'^(Table [12]\.|表[12]\s)',line):
            p.paragraph_format.keep_with_next=True
        if in_references:
            p.paragraph_format.line_spacing=1.05
            p.paragraph_format.space_after=Pt(4)
            p.paragraph_format.keep_together=True
            for r in p.runs:r.font.size=Pt(10 if kind=='supplement' else 10.5)
    figures=sorted(f for f in (OUT/'figures').glob('*.png') if f.name.startswith('FigureS') == (kind=='supplement'))
    expected=6 if kind=='supplement' and (ROOT/'repro_v3_cellchat_narrative.json').exists() else 5
    assert len(figures)==expected,len(figures)
    for f in figures:
        if kind!='supplement':doc.add_page_break()
        label=f.stem.split('_')[0].replace('FigureS','Supplementary Figure S').replace('Figure','Figure ')
        heading=doc.add_paragraph(label,style='Heading 1')
        if kind=='supplement':heading.paragraph_format.page_break_before=True
        p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
        with Image.open(f) as im:w,h=im.size
        # Preserve aspect ratio while leaving room for the figure heading.
        width=min(7.2,8.65*w/h)
        p.add_run().add_picture(str(f),width=Inches(width))
    doc.core_properties.title=raw.splitlines()[0][2:]
    doc.core_properties.author='朱四欢' if lang=='zh' else 'Sihuan Zhu'
    doc.core_properties.subject='Cross-cohort computational reproducibility study of public mouse ischemia data'
    doc.core_properties.keywords='cerebral ischemia; ligand receptor; LIANA; reproducibility; DecontX'
    path=OUT/f'{kind}_{lang}.docx';doc.save(path)
    return {'kind':kind,'language':lang,'docx':str(path),'references':len(order),'figures':len(figures),'tables':len(doc.tables)}

if __name__=='__main__':
    reports=[build(lang, kind) for lang in ['zh','en'] for kind in ['manuscript','supplement']]
    (OUT/'document_build_audit.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(reports,ensure_ascii=False))
