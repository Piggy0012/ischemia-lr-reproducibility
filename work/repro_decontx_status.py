"""Small-table inventory and convergence status; does not load expression matrices."""
from pathlib import Path
import json,sys
import pandas as pd
from repro_decontx import ROOT,OUT,summarize

def status(require_complete=False):
    summarize()
    expected=pd.read_csv(ROOT/'revision_sample_manifest.tsv',sep='\t')
    records=[]
    for row in expected.itertuples(index=False):
        p=OUT/row.dataset/row.sample/'audit.json'
        if not p.exists():continue
        a=json.loads(p.read_text())
        assert a['requested_parameters']['seed']==20260911
        assert a['requested_parameters']['background'] is None
        assert a['requested_parameters']['delta']==[10,10]
        assert a['requested_parameters']['estimateDelta'] is True
        assert a['python_verification']=={'barcode_order':True,'corrected_totals':True,'finite_nonnegative':True}
        assert a['decontX']=='1.10.0' and a['R']=='4.6.1'
        records.append({'dataset':row.dataset,'sample':row.sample,
                        'condition':'Sham' if ('sham' in row.source_prefix.lower() or '_WTC' in row.source_prefix) else 'MCAO',
                        'last_logged_iteration':a['last_logged_iteration'],
                        'last_logged_max_divergence':a['last_logged_max_divergence'],
                        'convergence_threshold_reached':a['convergence_threshold_reached'],
                        'n_zero_corrected_total':a['n_zero_corrected_total'],
                        'fractional_values':a['fractional_values']})
    if records:
        d=pd.read_csv(OUT/'sample_contamination_summary.tsv',sep='\t')
        s=pd.DataFrame(records)
        d=d.merge(s,on=['dataset','sample'],validate='one_to_one')
        d.to_csv(OUT/'sample_contamination_summary_with_convergence.tsv',sep='\t',index=False)
        report={'expected_samples':len(expected),'completed_samples':len(d),
                'completed_cells':int(d.n_cells.sum()),
                'all_completed_converged':bool(d.convergence_threshold_reached.all()),
                'n_zero_corrected_total':int(d.n_zero_corrected_total.sum()),
                'all_matrices_fractional':bool(d.fractional_values.all()),
                'parameters_versions_and_matrix_audits_checked':True,
                'samples':d.astype(object).where(pd.notna(d),None).to_dict(orient='records')}
        (OUT/'completion_status.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(d[['sample','n_cells','median_estimated_contamination','removed_total_fraction','last_logged_iteration','convergence_threshold_reached']].to_string(index=False))
    if require_complete:assert len(records)==len(expected),'Not all expected sample audits are complete'

if __name__=='__main__':status('--require-complete' in sys.argv)
