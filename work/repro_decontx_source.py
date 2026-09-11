"""Select corrected matrices only by the documented convergence rule.

Default outputs are immutable. Only a default nonconvergent fit is replaced
by the separately retained uniform maxIter=2000 numerical follow-up. This
choice never depends on ligand-receptor outputs.
"""
from pathlib import Path
import json, hashlib
ROOT=Path(__file__).resolve().parent

def resolve(acc,sample):
    base=ROOT/'repro_decontx'/acc/sample
    assert (base/'audit.json').exists(),f'Default fit incomplete: {sample}'
    original=json.loads((base/'audit.json').read_text(encoding='utf-8'))
    extended=not bool(original['convergence_threshold_reached'])
    selected=ROOT/'repro_decontx_maxiter2000'/acc/sample if extended else base
    assert (selected/'audit.json').exists(),f'Numerical extension pending: {sample}'
    d=json.loads((selected/'audit.json').read_text(encoding='utf-8'))
    assert d['sample']==sample and d['dataset']==acc
    if extended:
        assert d['requested_parameters']['maxIter']==2000
        assert d['input']['input_sha256']==original['input']['input_sha256']
        for k in ['z_source','background','batch','seed','delta','estimateDelta','convergence']:
            assert d['requested_parameters'][k]==original['requested_parameters'][k],k
    return selected,{'selected_relative_directory':selected.relative_to(ROOT).as_posix(),
        'default_converged':bool(original['convergence_threshold_reached']),
        'numerical_extension_selected':extended,
        'selected_converged':bool(d['convergence_threshold_reached']),
        'selected_max_iter':d['requested_parameters']['maxIter'],
        'selected_last_divergence':d['last_logged_max_divergence'],
        'selected_audit_sha256':hashlib.sha256((selected/'audit.json').read_bytes()).hexdigest(),
        'default_audit_sha256':hashlib.sha256((base/'audit.json').read_bytes()).hexdigest()}
