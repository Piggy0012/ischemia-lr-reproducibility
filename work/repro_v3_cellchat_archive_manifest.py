"""Freeze the CellChat extension overlay; never include runtimes or binary packages."""
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json
ROOT=Path(__file__).resolve().parent
PROJECT=ROOT.parent
OUT=PROJECT/'outputs/reproducibility_v3_cellchat'
OUT.mkdir(parents=True,exist_ok=True)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8'))

summary=read(ROOT/'repro_v3_cellchat/tables/summary_audit.json')
assert summary['status']=='complete' and summary['n_samples']==11 and summary['n_target_keys']==32
assert read(ROOT/'repro_v3_cellchat_editorial_audit.json')['status']=='complete'
telemetry=list((ROOT/'repro_v3_cellchat/controlled').glob('*/*/run_telemetry.json'))
assert len(telemetry)==11 and all(read(p)['status']=='complete' for p in telemetry)
for p in telemetry:
    loss=read(p.parent/'lossless_export_audit.json')
    assert loss['exact_probability_roundtrip'] and loss['exact_pvalue_roundtrip']
for label in ['scientific_input_file_sha256','output_file_sha256']:
    for rel,expected in summary[label].items():assert sha(ROOT/rel)==expected,(label,rel)
review=read(OUT/'figure_source_data/figureS6_visual_review.json')
assert review['passed'] and review['preview_sha256']==sha(OUT/'figures/FigureS6_cellchat_controlled_extension.png')

files=set(p for p in ROOT.glob('repro_v3_cellchat*') if p.is_file())
files.add(ROOT/'repro_v3_independent_framework_plan.md')
for folder in ['repro_v3_cellchat','repro_v3_cellchat_probe','repro_v3_cellchat_environment']:
    for p in (ROOT/folder).rglob('*'):
        if p.is_file() and not any(x in {'rlibrary','downloads','__pycache__'} for x in p.relative_to(ROOT/folder).parts):files.add(p)
for folder in ['figures','figure_source_data']:
    files.update(p for p in (OUT/folder).rglob('*') if p.is_file())
# These are already part of the baseline project; include them explicitly so
# the new score summary also has its expression inputs bound in this overlay.
for rel in summary['scientific_input_file_sha256']:
    p=ROOT/rel
    if p.parts[len(ROOT.parts)]=='processed':files.add(p)

rows=[]
for p in sorted(files):
    rel=p.relative_to(PROJECT).as_posix()
    assert p.suffix.lower() not in {'.exe','.dll','.zip','.npz','.pyc'},rel
    assert p.stat().st_size<100*1024**2,rel
    if '/repro_v3_cellchat_probe/' in rel:
        role='upstream_source_resource_and_provenance'
        license_note='CellChat-derived source/database retains upstream GPL-3; see THIRD_PARTY_NOTICE.md and LICENSE.CellChat; analysis provenance records remain separately attributable'
    elif '/controlled/' in rel:role='controlled_model_input_output_or_execution_audit';license_note='Analysis artifact; native CellChat source and database licensing is preserved separately'
    elif '/repro_v3_cellchat_environment/' in rel:role='effective_environment_lock_and_restore_audit';license_note='Package metadata retains its package attribution; no package binaries included'
    elif '/figure' in rel:role='supplementary_figure_or_source_data_and_QA';license_note='Original analysis artifact'
    elif '/processed/' in rel:role='baseline_expression_summary_input';license_note='Original analysis artifact from cited GEO data'
    elif p.suffix in {'.py','.R'}:role='analysis_or_audit_code';license_note='Original project analysis code unless explicitly identified as an upstream copy'
    else:role='analysis_table_narrative_or_provenance';license_note='Original analysis artifact with cited external source attribution where applicable'
    rows.append({'path':rel,'bytes':p.stat().st_size,'sha256':sha(p),'role':role,'license_note':license_note})

excluded=['Unpacked R/Python runtimes and installed library trees','work/repro_v3_cellchat/downloads (42 package ZIPs)',
 'work/repro_v3_release_assets/cellchat_environment_windows_R461.zip (release/Zenodo asset only)',
 'work/repro_v3_release_assets/CellChat_official_source_*.tar.gz (full upstream source release asset)',
 'Original large gene-count caches and raw matrices; normalized 865-gene sparse model inputs are included',
 'This manifest and its TSV companion are not self-hashed']
manifest={'status':'frozen_complete_extension_overlay','created_at_utc':datetime.now(timezone.utc).isoformat(),
  'scope':'Native CellChat controlled shared-resource arm only; raw/primary original eleven libraries; no additional inferential disease P values',
  'file_count':len(rows),'total_bytes':sum(r['bytes'] for r in rows),'largest_file':max(rows,key=lambda r:r['bytes']),
  'relative_path_base':'Project repository root, preserving work/ and outputs/ paths',
  'runtime_requirement':'Use separately restored matching scientific Python and R 4.6.1; environment/restore manifests are included, binary asset belongs only in release storage',
  'model_execution':{'libraries_complete':11,'maximum_measured_process_tree_RSS_bytes':max(read(p)['peak_process_tree_rss_bytes'] for p in telemetry),
     'total_individual_model_wall_seconds':sum(read(p)['elapsed_wall_seconds'] for p in telemetry)},
  'summary_audit_sha256':sha(ROOT/'repro_v3_cellchat/tables/summary_audit.json'),
  'environment_manifest_sha256':sha(ROOT/'repro_v3_cellchat_environment/environment_manifest.json'),
  'all_native_RDS_and_full_network_zero_scores_internal_P_preserved':True,
  'normalized_sparse_inputs_included':True,'third_party_license_preserved':True,
  'separate_release_assets':[{'path':p.relative_to(PROJECT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)}
      for p in [ROOT/'repro_v3_release_assets/cellchat_environment_windows_R461.zip',
                PROJECT/read(ROOT/'repro_v3_cellchat_probe/corresponding_source_asset.json')['local_path']]],
  'exclusions':excluded,'files':rows}
dest=OUT/'cellchat_archive_manifest.json'
dest.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
with (OUT/'cellchat_archive_manifest.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({k:manifest[k] for k in ['status','file_count','total_bytes','largest_file','model_execution']},indent=2))
