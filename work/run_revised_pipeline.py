"""Rebuild revision outputs from archived sufficient statistics or submitted matrices."""
from pathlib import Path
import argparse,subprocess,sys,os,shutil,tarfile
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--recompute-statistics',action='store_true',help='refit revised pseudobulk and custom scores from included aggregates')
p.add_argument('--from-submitted-matrices',action='store_true',help='download public matrices if absent; force new annotations, aggregation, LIANA and revised statistics')
p.add_argument('--document-python',default=sys.executable,help='Python with python-docx; use the bundled runtime in Codex')
args=p.parse_args()
def run(script,*extra,python=None):
    subprocess.run([python or sys.executable,'-X','utf8',str(ROOT/script),*extra],check=True,cwd=ROOT.parent)
if args.from_submitted_matrices:
    os.environ['ISCHEMIA_REVISION_FORCE']='1'
    import pandas as pd
    manifest=pd.read_csv(ROOT/'revision_sample_manifest.tsv',sep='\t')
    for acc in ['GSE174574','GSE245386']:
        expected=manifest[manifest.dataset==acc]
        if not all((ROOT/f'data/{acc}/raw/{prefix}_matrix.mtx.gz').exists() for prefix in expected.source_prefix):
            if acc=='GSE174574':
                run('download_geo.py',acc,'--matrix')
                raw=ROOT/f'data/{acc}/raw';raw.mkdir(parents=True,exist_ok=True)
                with tarfile.open(ROOT/f'data/{acc}/{acc}_RAW.tar') as archive:archive.extractall(raw,filter='data')
            else:run('download_validation.py')
        # The frozen original rule labels and feature tables are included. If
        # absent, recreate them from the same submitted count matrices.
        if not all((ROOT/f'processed/{acc}/{sample}_cells.tsv.gz').exists() for sample in expected['sample']):run('stream_annotate.py',acc)
        run('revision_identity.py',acc);run('revision_aggregate.py',acc);run('revision_liana.py',acc)
    run('revision_gene_audit.py','markers')
if args.recompute_statistics or args.from_submitted_matrices:
    os.environ['ISCHEMIA_FORCE_RECOMPUTE']='1';run('revision_statistics.py')
run('revision_summary.py');run('revision_gene_audit.py','bbb')
run('revision_identity_figures.py');run('revision_robustness_figures.py','--export');run('revision_barrier_figure.py')
out=ROOT.parent/'outputs/revised';(out/'figures').mkdir(parents=True,exist_ok=True)
for f in (ROOT.parent/'outputs/figures').glob('*'):
    if f.is_file():shutil.copy2(f,out/'figures'/f.name)
run('build_revised_manuscript.py',python=args.document_python)
run('validate_revision.py')
print('Revision rebuilt in',out)
