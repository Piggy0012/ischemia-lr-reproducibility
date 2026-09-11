"""Download and stream-QC the independently deposited GSE332910 RNA libraries.

Library replication is established; individual-mouse versus pool provenance is
not established by the accessible source metadata. Do not label these as six
independent animals. Uses the ORIGINAL stream_annotate lineage rule unchanged.
The deposited PISA matrix is gene-major, unlike the original 10x matrices.
"""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, gzip, hashlib, json, shutil, tarfile, time, urllib.request
import pandas as pd
import stream_annotate

BASE = Path(__file__).resolve().parent
ROOT = BASE / 'repro_third_cohort'
ACC = 'GSE332910'
SAMPLES = [
    ('GSM9755192', 'MCAO', 1, 'snRNA-seq_24h_1.tar.gz', 137766625),
    ('GSM9755194', 'MCAO', 2, 'snRNA-seq_24h_2.tar.gz', 101499743),
    ('GSM9755196', 'MCAO', 3, 'snRNA-seq_24h_3.tar.gz', 93333915),
    ('GSM9755210', 'Sham', 1, 'snRNA-seq_sham_1.tar.gz', 88263041),
    ('GSM9755212', 'Sham', 2, 'snRNA-seq_sham_2.tar.gz', 101326327),
    ('GSM9755214', 'Sham', 3, 'snRNA-seq_sham_3.tar.gz', 39760025),
]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def download(row):
    gsm, condition, replicate, suffix, size = row
    name = gsm + '_' + suffix
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9755nnn/{gsm}/suppl/{name}'
    dest = ROOT / 'archives' / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists() or dest.stat().st_size != size:
        partial = dest.with_suffix(dest.suffix + '.partial')
        for attempt in range(5):
            try:
                # Restart partial downloads: GEO does not always honour Range.
                with urllib.request.urlopen(url, timeout=60) as response, partial.open('wb') as f:
                    total = 0
                    while True:
                        b = response.read(1024 * 1024)
                        if not b:
                            break
                        f.write(b)
                        total += len(b)
                if total != size:
                    raise ValueError(f'{gsm}: downloaded {total}, expected {size}')
                partial.replace(dest)
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(2 + attempt)
    raw = ROOT / 'data' / ACC / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    prefix = f'{gsm}_{condition}{replicate}'
    extracted = []
    for member in tarfile.open(dest, mode='r:gz'):
        name = Path(member.name).name
        if not member.isfile() or name.startswith('._'):
            continue
        kind = next((x for x in ('matrix.mtx.gz', 'features.tsv.gz', 'barcodes.tsv.gz') if name.endswith(x)), None)
        if kind is None:
            continue
        dst = raw / (prefix + '_' + kind)
        if kind == 'features.tsv.gz':
            # The source has one column of submitted mouse gene symbols, NOT Ensembl IDs.
            with tarfile.open(dest, mode='r:gz') as archive:
                symbols = gzip.decompress(archive.extractfile(member.name).read()).decode().splitlines()
            assert symbols and all('\t' not in s for s in symbols)
            assert len(set(symbols)) == len(symbols), 'Duplicate submitted symbols require explicit aggregation.'
            with gzip.open(dst, 'wt', encoding='utf-8', newline='') as f:
                for symbol in symbols:
                    f.write(symbol + '\t' + symbol + '\n')
            feature_report = {'n_genes': len(symbols), 'n_unique_symbols': len(set(symbols)),
                              'first_symbols': symbols[:5], 'mapping': 'Source one-column gene symbols copied as both submitted identifier and symbol; no orthology mapping or invented Ensembl IDs.'}
        elif not dst.exists() or dst.stat().st_size != member.size:
            with tarfile.open(dest, mode='r:gz') as archive, dst.open('wb') as f:
                shutil.copyfileobj(archive.extractfile(member.name), f, 1024 * 1024)
        extracted.append({'source_member': member.name, 'destination_relative': str(dst.relative_to(ROOT)), 'size_bytes': dst.stat().st_size, 'sha256': digest(dst)})
    assert len(extracted) == 3, (gsm, extracted)
    with gzip.open(raw / (prefix + '_matrix.mtx.gz'), 'rt') as f:
        header = f.readline().strip()
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        shape = list(map(int, line.split()))
        first_entries = [f.readline().strip() for _ in range(3)]
    record = {'accession': ACC, 'sample': gsm, 'condition': condition, 'library_replicate': replicate,
              'n_independent_animals': None, 'pooling_status': 'individual_mouse_or_pool_mapping_not_explicitly_reported',
              'author_reported_biological_replicates_per_condition': 3,
              'biological_replicate_source': 'Supplementary Table S1, Sheet1 rows 2-4, ADVS-9999-e77547-s001.xlsx; https://pmc-oa-opendata.s3.amazonaws.com/PMC13542395.1/ADVS-9999-e77547-s001.xlsx',
              'time': '24 h after 60-minute MCAO reperfusion' if condition == 'MCAO' else 'Sham surgery',
              'tissue': 'ipsilateral striatum', 'assay': 'DNBelab single-nucleus RNA sequencing; MGI DNBSEQ-Tx',
              'source_url': url, 'geo_url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}',
              'publication_doi': '10.1002/advs.77547', 'archive_sha256': digest(dest), 'archive_bytes': size,
              'matrix_header': header, 'matrix_shape_genes_cells_nnz': shape, 'first_entries': first_entries,
              'feature_mapping': feature_report, 'extracted_files': extracted,
              'metadata_notes': ['The series-level unimmunized/immunized wording conflicts with the stroke design; group assignment uses explicit sample title/treatment and source article.',
                                 'ATAC and RNA records of the same named replicate are not independent biological samples.',
                                 'Counts are deposited processed integer UMI matrices, not an assertion of empty-droplet availability.']}
    (ROOT / (gsm + '_source_audit.json')).write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
    print('READY', gsm, shape, flush=True)
    return record

def qc():
    stream_annotate.ROOT = ROOT
    stream_annotate.process(ACC)
    out = ROOT / 'processed' / ACC
    q = pd.read_csv(out / 'qc_summary.tsv', sep='\t')
    records = [json.loads((ROOT / (row[0] + '_source_audit.json')).read_text(encoding='utf-8')) for row in SAMPLES]
    r = pd.DataFrame([{'sample': x['sample'], 'condition': x['condition'], 'library_replicate': x['library_replicate'], 'n_independent_animals': None, 'n_raw_genes': x['matrix_shape_genes_cells_nnz'][0], 'n_raw_cells': x['matrix_shape_genes_cells_nnz'][1], 'n_raw_nonzero': x['matrix_shape_genes_cells_nnz'][2]} for x in records])
    q['n_qc_unassigned'] = q.n_qc - q.n_assigned
    q['qc_unassigned_fraction'] = q.n_qc_unassigned / q.n_qc
    q['author_reported_biological_replicates_per_condition'] = 3
    q.merge(r, on='sample', validate='1:1').to_csv(ROOT / 'library_qc_summary.tsv', sep='\t', index=False)
    (ROOT / 'download_manifest.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
    config = {'qc': '200-6000 detected genes inclusive; >=500 UMI; <=20% mitochondrial UMI',
              'annotation': 'Unmodified work/stream_annotate.py marker panels and scoring, calculated per library.',
              'source_script_sha256': digest(BASE / 'stream_annotate.py'),
              'status': 'Exploratory third cohort. Author Table S1 reports three biological replicates and three RNA libraries per group; individual mouse/pool mapping unresolved.',
              'comparison_notes': 'snRNA-seq, striatum and BGI chemistry differ from whole-hemisphere scRNA-seq cohorts. This is not a replacement for same-assay validation.',
              'python_version': __import__('sys').version, 'pandas_version': pd.__version__, 'numpy_version': __import__('numpy').__version__}
    (ROOT / 'qc_method.json').write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding='utf-8')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--download-only', action='store_true')
    ap.add_argument('--qc-only', action='store_true')
    args = ap.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    if not args.qc_only:
        records = [download(row) for row in SAMPLES]
        (ROOT / 'download_manifest.json').write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding='utf-8')
    if not args.download_only:
        qc()
