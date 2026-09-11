"""Re-use the original, sample-level statistical routines on revised selections."""
from pathlib import Path
import sys,os
import pandas as pd
import run_pseudobulk_de as de
import candidate_communication as communication
ROOT=Path(__file__).resolve().parent
de.ROOT=ROOT/'revision_analysis';communication.ROOT=de.ROOT
for acc in sys.argv[1:] or ['GSE174574','GSE245386']:
    counts=pd.read_csv(de.ROOT/f'processed/{acc}/sample_cell_counts.tsv',sep='\t')
    for config in ['singlet','reference_singlet','reference_only']:
        targets=counts[(counts.config==config)&counts.cell_type.isin(['Astrocyte','Endothelial'])]
        assert len(targets)==(12 if acc=='GSE174574' else 10) and targets.n_cells.ge(30).all(), 'Insufficient target cells: define a matched DE/communication sample universe before inference'
        de.process(acc,[config])
        if os.environ.get('ISCHEMIA_FORCE_RECOMPUTE')=='1' or not (de.ROOT/f'results/{acc}/{config}__communication.tsv.gz').exists():communication.process(acc,config)
