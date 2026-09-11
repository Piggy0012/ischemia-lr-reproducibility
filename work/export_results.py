"""Copy the final numerical results to the user-facing output tree."""
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'outputs'
OUT.mkdir(exist_ok=True)
shutil.copytree(ROOT/'results',OUT/'tables',dirs_exist_ok=True)
sorted_out=OUT/'tables/sorted_validation'
sorted_out.mkdir(exist_ok=True)
for name in ['diff_expression.tsv.gz','sample_manifest.tsv','audit.json','identity_markers_expression.tsv']:
    source=ROOT/'sorted_rna'/name
    if source.exists():shutil.copy2(source,sorted_out/name)
shutil.copy2(ROOT/'dataset_audit.json',OUT/'dataset_audit.json')
shutil.copy2(ROOT/'delivery_readme.md',OUT/'README.md')
print('Exported result tables and reading guide.')
