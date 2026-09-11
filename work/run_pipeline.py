"""Reproduce the deposited-matrix analysis or reuse bundled sufficient statistics."""
from pathlib import Path
import argparse,subprocess,sys,tarfile,json,gzip,os,shutil,datetime
ROOT=Path(__file__).resolve().parent
def run(script,*args):subprocess.run([sys.executable,'-X','utf8',str(ROOT/script),*args],check=True,cwd=ROOT.parent)
p=argparse.ArgumentParser();p.add_argument('--from-raw',action='store_true',help='download submitted matrices and recompute QC/annotation in a new isolated directory (network and approximately 4 GB disk required)');p.add_argument('--use-cache',action='store_true',help='reuse existing DE/communication tables; rebuild downstream summaries and figures');p.add_argument('--raw-worker',action='store_true',help=argparse.SUPPRESS);args=p.parse_args()
if args.from_raw and not args.raw_worker:
 fresh=ROOT.parent/('fresh_run_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))/'work';fresh.mkdir(parents=True,exist_ok=False)
 for f in ROOT.glob('*.py'):shutil.copy2(f,fresh/f.name)
 for f in ROOT.glob('*.md'):shutil.copy2(f,fresh/f.name)
 for f in ROOT.glob('*.json'):shutil.copy2(f,fresh/f.name)
 for folder in ['literature','metadata_audit','spatial_geometry','sorted_rna']:shutil.copytree(ROOT/folder,fresh/folder)
 subprocess.run([sys.executable,'-X','utf8',str(fresh/'run_pipeline.py'),'--from-raw','--raw-worker'],check=True,cwd=fresh.parent);sys.exit(0)
if not args.use_cache:os.environ['ISCHEMIA_FORCE_RECOMPUTE']='1'
else:os.environ.pop('ISCHEMIA_FORCE_RECOMPUTE',None)
(ROOT/'results').mkdir(exist_ok=True)
if args.from_raw:
 run('download_geo.py','GSE174574','--matrix')
 raw=ROOT/'data/GSE174574/raw';raw.mkdir(parents=True,exist_ok=True)
 with tarfile.open(ROOT/'data/GSE174574/GSE174574_RAW.tar') as t:t.extractall(raw,filter='data')
 run('download_validation.py')
 for acc in ['GSE174574','GSE245386']:
  run('stream_annotate.py',acc);run('aggregate_targets.py',acc);run('contamination_sensitivity.py',acc)
 run('identity_audit.py');run('download_spatial.py');run('spatial_analysis.py')
run('sorted_rna/analyze_sorted.py')
for acc in ['GSE174574','GSE245386']:run('run_pseudobulk_de.py',acc,'primary','strict_identity','mt10','low_myeloid')
run('run_all_communication.py');run('cross_cohort.py');run('integrated_tables.py');run('make_figures.py');run('build_manuscript.py');run('validate_results.py')
run('export_results.py')
print('Completed. Outputs are in',ROOT.parent/'outputs')
