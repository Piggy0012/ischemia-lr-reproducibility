"""Capture effective package versions, respecting import-path precedence."""
from pathlib import Path
import importlib.metadata as metadata
import sys,platform,json,re
ROOT=Path(__file__).resolve().parent
names={}
for d in metadata.distributions():
    name=d.metadata.get('Name')
    if name:
        canonical=re.sub(r'[-_.]+','-',name).lower()
        names.setdefault(canonical,name)
# metadata.version follows effective path precedence. A dictionary populated
# from all distributions can otherwise silently overwrite a venv version
# with a later, shadowed system-site distribution.
packages={name:metadata.version(name) for name in sorted(names)}
import pandas,numpy,scipy
effective={'pandas':pandas.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__}
assert all(packages[k]==v for k,v in effective.items())
report={'python':platform.python_version(),'platform':platform.platform(),'packages':packages,
    'import_verified_versions':effective,'environment':'revision_env with system-site access; effective versions resolve first on sys.path',
    'metadata_resolution':'importlib.metadata.version per canonical name; do not overwrite earlier distributions',
    'lockfile':'requirements_reproducibility_lock.txt',
    'excluded_from_portable_lock':[{'package':'artifact-tool-v2','reason':'Host-only local artifact helper is not used by the scientific analysis; its local direct URL is not portable.'}],
    'provenance_correction':'Earlier metadata inventory listed a shadowed pandas 3.0.1 from system sites. Analysis imports and the freeze lock use pandas 2.3.3; this capture records effective versions.'}
(ROOT/'repro_environment.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({'python':report['python'],'effective_imports':effective,'packages':len(packages)}))
