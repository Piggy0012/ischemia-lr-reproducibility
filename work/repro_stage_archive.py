"""Prepare a self-contained code/results release, keeping raw downloads separate.

Does not publish or invent a DOI. Run after the v2 manuscript and audits are final.
"""
from pathlib import Path
import shutil,json,hashlib,sys
ROOT=Path(__file__).resolve().parent
STAGE=ROOT/'repro_delivery/analysis'
OLD=ROOT/'revised_delivery/analysis'
OUT=ROOT.parent/'outputs/reproducibility_v2'
assert (OLD/'work/run_revised_pipeline.py').exists()
STAGE.mkdir(parents=True,exist_ok=True)

def copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)

def tree(src,dst,exclude=()):
    if src.exists():shutil.copytree(src,dst,dirs_exist_ok=True,
        ignore=shutil.ignore_patterns('__pycache__','*.pyc',*exclude))

# The previously audited stage contains code, cell summaries and all tests,
# not installed runtimes or multi-gigabyte submitted matrices.
tree(OLD/'work',STAGE/'work')
# The current descriptive audit hashes this exact compressed source. The older
# stage has identical decoded values but a different gzip header; preserve the
# current audited bytes instead of silently inheriting the older container.
copy(ROOT/'sorted_rna/astro_endothelial_normalized.tsv.gz',
     STAGE/'work/sorted_rna/astro_endothelial_normalized.tsv.gz')
for pattern in ['repro_*.py','repro_*.md','repro_*.json','repro_*.R','requirements_revision.txt','requirements_reproducibility_lock.txt','software_versions_revision.json']:
    for f in ROOT.glob(pattern):copy(f,STAGE/'work'/f.name)
for folder in ['repro_literature','repro_decontx_analysis','repro_liana_rank_diagnostic','repro_third_liana']:
    tree(ROOT/folder,STAGE/'work'/folder)
for f in (ROOT/'repro_cohort_metadata').glob('*'):
    if f.suffix in ['.json','.soft','.xml','.xlsx'] or f.name.endswith('_paper.txt'):
        copy(f,STAGE/'work/repro_cohort_metadata'/f.name)
third=ROOT/'repro_third_cohort'
for folder in ['processed','results']:tree(third/folder,STAGE/'work/repro_third_cohort'/folder)
for f in third.glob('*'):
    if f.is_file() and f.suffix in ['.json','.tsv','.md']:copy(f,STAGE/'work/repro_third_cohort'/f.name)
for sample in (third/'context_cache').glob('GSM*'):
    for f in sample.glob('*'):
        if f.is_file() and (f.suffix in ['.json','.tsv'] or f.name.endswith('.tsv.gz')):
            copy(f,STAGE/'work/repro_third_cohort/context_cache'/sample.name/f.name)
# Preserve model and software audits without installing R or copying dependencies.
de=ROOT/'repro_decontx'
tree(de/'source_snapshots',STAGE/'work/repro_decontx/source_snapshots')
for pattern in ['*.json','*.tsv','*.md']:
    for f in de.glob(pattern):copy(f,STAGE/'work/repro_decontx'/f.name)
for f in (de/'installers').glob('*'):
    if f.suffix in ['.json','.R','.py','.txt'] or f.name.endswith(('_LICENSE','_DESCRIPTION')):copy(f,STAGE/'work/repro_decontx/installers'/f.name)
for source in ['repro_decontx','repro_decontx_maxiter2000']:
    for f in (ROOT/source).glob('*'):
        if f.is_file() and f.suffix in ['.json','.tsv','.md','.R','.py','.txt']:
            copy(f,STAGE/'work'/source/f.name)
    for acc in ['GSE174574','GSE245386']:
        for sample in (ROOT/source/acc).glob('GSM*'):
            for f in sample.glob('*'):
                if f.is_file() and (f.suffix in ['.json','.R','.py','.log','.txt','.tsv','.rds'] or f.name.endswith('.tsv.gz')):
                    copy(f,STAGE/'work'/source/acc/sample.name/f.name)
tree(OUT,STAGE/'outputs/reproducibility_v2',exclude=('*.zip','release_integrity.json'))
# Historical outputs remain labeled, but old manuscript drafts are not the entry point.
tree(OLD/'outputs/tables',STAGE/'outputs/tables')
tree(OLD/'outputs/figure_source_data',STAGE/'outputs/figure_source_data')
tree(OLD/'outputs/revised/tables',STAGE/'outputs/revised/tables')
tree(OLD/'outputs/revised/figure_source_data',STAGE/'outputs/revised/figure_source_data')

author=json.loads((ROOT/'repro_author.json').read_text(encoding='utf-8'))
title='Expression-level concordance and implementation sensitivity of ligand–receptor priority rankings in mouse cerebral ischemia'
record=OUT/'archive_record.json'
published=json.loads(record.read_text(encoding='utf-8')) if record.exists() else {}
doi=published.get('doi')
metadata={'title':title,'upload_type':'software','creators':[{'name':author.get('archive_name',author['name']),'orcid':author['orcid'],'affiliation':author['affiliation']}],
    'description':'Code, environment records, pseudobulk counts, complete tests and label permutations, full ligand–receptor network scores, implementation diagnostics, and ambient-RNA sensitivity outputs supporting an exploratory cross-cohort reproducibility assessment. Native and diagnostic aggregation and raw versus estimated decontaminated counts are distinguished. Original GEO downloads remain traceable by accession and checksum.',
    'access_right':'open','license':'MIT','version':'2.0.0',
    'keywords':['cerebral ischemia','cell-cell communication','reproducibility','ligand receptor','LIANA','DecontX']}
if doi:metadata['doi']=doi
(STAGE/'.zenodo.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
cff='cff-version: 1.2.0\nmessage: "Please cite this versioned analysis archive and the original GEO studies."\ntype: software\n'
cff+='title: '+json.dumps(title,ensure_ascii=False)+'\nversion: 2.0.0\nauthors:\n'
cff+='  - family-names: '+json.dumps(author.get('family_name_en',author['family_name']),ensure_ascii=False)+'\n    given-names: '+json.dumps(author.get('given_name_en',author['given_name']),ensure_ascii=False)+'\n    affiliation: '+json.dumps(author['affiliation'],ensure_ascii=False)+'\n    orcid: "https://orcid.org/'+author['orcid']+'"\n'
if doi:cff+='doi: "'+doi+'"\n'
(STAGE/'CITATION.cff').write_text(cff,encoding='utf-8')
(STAGE/'LICENSE').write_text('''MIT License

Copyright (c) 2026 朱四欢

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This license applies to the original analysis code in this release. Original
datasets, interaction resources, pretrained models and copied upstream source
retain their original terms and attribution; this release does not relicense
third-party materials. Original derived tables, manuscript and figures are
shared under Creative Commons Attribution 4.0 International (CC BY 4.0).
''',encoding='utf-8')
(STAGE/'README.md').write_text('''# Cross-cohort ligand–receptor reproducibility analysis

The current manuscript and revision response are in `outputs/reproducibility_v2/`.
The author supplied for this archive is 朱四欢, ORCID 0009-0005-9705-0585.

Start with `outputs/reproducibility_v2/README.md`. Tables with `null_`, `rank_`,
`ambient_`, `third_`, and `sorted_` prefixes identify distinct analyses. Original
and earlier revision tests remain in `work/results`, `work/revision_analysis`,
and `work/sorted_rna`; their existence does not make them confirmatory endpoints.

The code/results archive includes original integer UMI pseudobulk counts and
all statistical results. Fractional DecontX estimates are labeled separately.
Large corrected single-cell matrices, if distributed, are in the separate
matrix archive. Original submitted matrices are obtained using archived GEO
URLs and checksums. No private account credentials or installed runtimes are
part of the release.

The packaged LIANA diagnostics do not modify the installed library and do not
claim to be an official patched release. Full per-library network scores are
retained so relative-ranking diagnostics can be rerun without a raw-cell fit.

Archive metadata alone is not proof of public deposition. The final public
record and DOI, when published, are recorded in `archive_record.json`.
''',encoding='utf-8')
print('Staged code/results at',STAGE)
