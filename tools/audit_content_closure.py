"""Check original payload preservation, named key data and repository-entry links."""
from pathlib import Path
import csv
import gzip
import json
import re
from update_file_index import sha, eligible

ROOT = Path(__file__).resolve().parent.parent


def main():
    original = list(csv.DictReader((ROOT / 'MANIFEST_SHA256.tsv').open(encoding='utf-8'), delimiter='\t'))
    replaced = {'README.md','CITATION.cff','.zenodo.json'}
    for row in original:
        path = ROOT / ('provenance/v2_archive/' + row['path'] if row['path'] in replaced else row['path'])
        assert path.exists() and path.stat().st_size == int(row['bytes']) and sha(path) == row['sha256'], row['path']
    env = ROOT / 'work/repro_v3_environment'
    lock = json.loads((env / 'environment_manifest.json').read_text(encoding='utf-8'))
    for row in lock['file_records']:
        p = env / row['path']
        assert p.stat().st_size == row['bytes'] and sha(p) == row['sha256'], row['path']
    key_paths = [
        'work/processed/GSE174574/marker_panel.json',
        'work/processed/GSE245386/marker_panel.json',
        'work/repro_third_cohort/processed/GSE332910/marker_panel.json',
        'outputs/reproducibility_v2/tables/third_rank_diagnostic_original_pair_on99.tsv',
        'outputs/reproducibility_v2/tables/third_rank_diagnostic_original_pair_on99_effects.tsv.gz',
        'outputs/reproducibility_v2/tables/third_rank_diagnostic_original_pair_on99_audit.json',
        'outputs/reproducibility_v2/tables/sorted_all_gene_descriptive_effects.tsv.gz',
    ]
    for p in key_paths:
        assert (ROOT / p).is_file(), p
    with gzip.open(ROOT / key_paths[-1], 'rt', encoding='utf-8', newline='') as f:
        n_effect_rows = sum(1 for _ in csv.reader(f, delimiter='\t')) - 1
    assert n_effect_rows == 49914
    links = []
    for name in ['README.md','RUNNING.md','ENVIRONMENT.md','AUTHOR_DECLARATIONS.md']:
        text = (ROOT / name).read_text(encoding='utf-8')
        for value in re.findall(r'\]\(([^)]+)\)', text):
            if '://' in value or value.startswith('#'):
                continue
            dest = ROOT / value.split('#')[0]
            assert dest.exists(), (name, value)
            links.append({'document':name,'target':value})
    assert '* -text' in (ROOT / '.gitattributes').read_text()
    assert not any(p.stat().st_size >= 100 * 2**20 for p in eligible())
    result = {'status':'passed','original_v2_payload_files_verified':len(original),
              'original_scientific_payload_bytes_preserved':True,
              'updated_repository_metadata_with_exact_original_backups':sorted(replaced),
              'environment_payload_files_verified':len(lock['file_records'])+1,
              'all_named_key_data_present':True,'sorted_all_gene_descriptive_rows':n_effect_rows,
              'key_data':[{'path':p,'bytes':(ROOT/p).stat().st_size,'sha256':sha(ROOT/p)} for p in key_paths],
              'repository_entry_local_links_verified':links,
              'git_text_conversion_disabled':True,'files_at_least_100_MiB':[],
              'v3_manuscript_synchronized':False,'independent_framework_runtime_synchronized':False,
              'public_remote_existence_asserted':False}
    (ROOT / 'provenance/content_closure_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in {'key_data','repository_entry_local_links_verified'}},indent=2))


if __name__=='__main__':
    main()
