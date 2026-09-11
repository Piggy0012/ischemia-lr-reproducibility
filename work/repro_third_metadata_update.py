"""Record the successfully retrieved author replication metadata and QC denominator."""
from pathlib import Path
import hashlib, json
import pandas as pd
ROOT = Path(__file__).resolve().parent
COHORT = ROOT/'repro_third_cohort'
META = ROOT/'repro_cohort_metadata'
TABLES = ROOT.parent/'outputs/reproducibility_v2/tables'
SOURCE = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC13542395.1/ADVS-9999-e77547-s001.xlsx'
LOCATOR = 'Supplementary Table S1, Sheet1 B2:D4 (Sham and 24 h columns)'
table = json.loads((META/'ADVS-9999-e77547-s001.json').read_text(encoding='utf-8'))['Sheet1']
assert table[1][1] == 'Sham' and table[1][3] == '24 h'
assert table[2][0] == 'Biological replicates' and table[3][0] == 'Number of libraries'
assert table[2][1] == table[2][3] == table[3][1] == table[3][3] == 3
records = json.loads((COHORT/'download_manifest.json').read_text(encoding='utf-8'))
rows, qc_rows = [], []
for r in records:
    r['author_reported_biological_replicates_per_condition'] = 3
    r['biological_replicate_source'] = LOCATOR+'; '+SOURCE
    r['pooling_status'] = 'individual_mouse_or_pool_mapping_not_explicitly_reported'
    (COHORT/(r['sample']+'_source_audit.json')).write_text(json.dumps(r, ensure_ascii=False, indent=2), encoding='utf-8')
    rows.append({k:r[k] for k in ['sample','condition','library_replicate','n_independent_animals',
        'author_reported_biological_replicates_per_condition','biological_replicate_source','pooling_status',
        'time','tissue','assay','publication_doi','geo_url']})
    q = pd.read_csv(COHORT/'processed/GSE332910'/(r['sample']+'_cells.tsv.gz'),sep='\t')
    assert q.loc[q.cell_type != 'Unassigned','qc_pass'].all()
    qc_rows.append({'sample':r['sample'],'n_qc_unassigned':int((q.qc_pass & (q.cell_type=='Unassigned')).sum()),
        'qc_unassigned_fraction':float((q.loc[q.qc_pass,'cell_type']=='Unassigned').mean())})
(COHORT/'download_manifest.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
pd.DataFrame(rows).to_csv(TABLES/'third_library_metadata.tsv',sep='\t',index=False)
q = pd.read_csv(COHORT/'library_qc_summary.tsv',sep='\t').drop(columns=['n_qc_unassigned','qc_unassigned_fraction'],errors='ignore')
q = q.merge(pd.DataFrame(qc_rows),on='sample',validate='1:1')
q['author_reported_biological_replicates_per_condition'] = 3
assert (q.n_qc_unassigned == q.n_qc-q.n_assigned).all()
q.to_csv(COHORT/'library_qc_summary.tsv',sep='\t',index=False)
report = {'checked_date':'2026-09-11','dataset':'GSE332910','source_url':SOURCE,'source_locator':LOCATOR,
    'source_sha256':hashlib.sha256((META/'ADVS-9999-e77547-s001.xlsx').read_bytes()).hexdigest(),
    'author_reported_biological_replicates':{'Sham':3,'24 h':3},'reported_number_libraries':{'Sham':3,'24 h':3},
    'individual_mouse_id_or_pool_size_map':None,
    'interpretation':'Author-reported biological replication is verified in Table S1; one mouse per library or independent pool membership is not individually mapped.',
    'successful_access':'Public PMC Cloud Service; original metadata MD5 verified for both supplementary files.',
    'retired_endpoint_attempts':[{'url':u,'HTTP_status':404} for u in [
        'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC13542395',
        'https://pmc.ncbi.nlm.nih.gov/utils/oa/oa.fcgi?id=PMC13542395']],
    'official_service_transition_documentation':'https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/',
    'n_qc':int(q.n_qc.sum()),'n_qc_unassigned':int(q.n_qc_unassigned.sum()),
    'qc_unassigned_fraction':float(q.n_qc_unassigned.sum()/q.n_qc.sum()),
    'all_assigned_cells_confirmed_qc_pass':True}
(META/'GSE332910_replication_provenance_update.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
p = ROOT/'repro_external_cohort_audit.json'
a = json.loads(p.read_text(encoding='utf-8'))
a['analysis_completion_date'] = '2026-09-11'
a['selection'] = 'GSE332910 added as exploratory acute third cohort, with six downloaded and QC-verified libraries and author-reported three biological replicates per group (Table S1); individual mouse/pool mapping is not provided.'
a['third_cohort_replication_update'] = report
coavailability_path = TABLES/'third_coavailability_audit.json'
coavailability = json.loads(coavailability_path.read_text(encoding='utf-8'))
coavailability['author_reported_biological_replicates_per_condition'] = 3
coavailability['biological_replicate_source'] = SOURCE
coavailability_path.write_text(json.dumps(coavailability,indent=2),encoding='utf-8')
a['third_cohort_analysis'] = coavailability
for label,name in [('third_liana_analysis','third_rank_diagnostic_audit.json'),
                   ('original_pair_on_same_99','third_rank_diagnostic_original_pair_on99_audit.json')]:
    extra=TABLES/name
    if extra.exists():
        a[label]=json.loads(extra.read_text(encoding='utf-8'))
for c in a['new_candidates']:
    if c['accession']=='GSE332910':
        c['decision']='included_exploratory_acute_third_cohort'
        c['reason']='Table S1 reports three biological replicates and RNA libraries per Sham/24h condition; original matrices downloaded, QC complete and both target types adequate. Individual mouse/pool map absent; assay/platform/region differ.'
p.write_text(json.dumps(a,ensure_ascii=False,indent=2),encoding='utf-8')
qm = COHORT/'qc_method.json'
d = json.loads(qm.read_text())
d['status']='Exploratory third cohort. Author Table S1 reports three biological replicates and three RNA libraries per group; individual mouse/pool mapping unresolved.'
qm.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
