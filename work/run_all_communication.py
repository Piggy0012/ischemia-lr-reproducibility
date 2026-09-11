from pathlib import Path
from candidate_communication import process
for acc in ['GSE174574','GSE245386']:
 for config in ['primary','strict_identity','mt10','low_myeloid']:
  if __import__('os').environ.get('ISCHEMIA_FORCE_RECOMPUTE') or not (Path(__file__).resolve().parent/f'results/{acc}/{config}__communication.tsv.gz').exists():process(acc,config)
