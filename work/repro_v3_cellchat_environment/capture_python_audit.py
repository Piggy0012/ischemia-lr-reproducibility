"""Capture optional independent-audit dependencies without changing either runtime.

Baseline locks remain immutable. The current complete pip lock is supplementary;
the small delta is sufficient for an already matching baseline environment.
"""
from pathlib import Path
from datetime import datetime,timezone
import ast,hashlib,importlib,importlib.metadata as md,json,os,platform,re,subprocess,sys
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
BASE=ROOT/'work/repro_v3_environment/requirements-python.lock.txt'
SCRIPT=ROOT/'work/repro_v3_cellchat_independent_review.py'

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def locks(text):
    result={}
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):continue
        assert '==' in line and ' @ ' not in line,line
        name,version=line.split('==');key=canonicalize_name(name)
        assert key not in result,key
        result[key]={'name':name,'version':version,'line':line}
    return result

def main():
    assert sys.version_info[:3]==(3,12,14)
    env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONUTF8']='1'
    freeze=subprocess.check_output([sys.executable,'-m','pip','freeze','--all'],text=True,encoding='utf-8',env=env)
    excluded=[];kept=[]
    for line in freeze.splitlines():
        if line.lower().startswith(('artifact-tool-v2','artifact_tool_v2')):
            excluded.append({'package':'artifact-tool-v2','reason':'Host artifact tooling with machine-local source; absent from independent-review import/dependency closure.'})
        elif line.strip():kept.append(line)
    current=locks('\n'.join(kept));baseline=locks(BASE.read_text(encoding='utf-8'))
    added=sorted(set(current)-set(baseline));removed=sorted(set(baseline)-set(current))
    changed=[k for k in sorted(set(current)&set(baseline)) if current[k]['version']!=baseline[k]['version']]
    assert not removed and not changed,(removed,changed)
    (HERE/'requirements-python-current.lock.txt').write_text('\n'.join(current[k]['line'] for k in sorted(current))+'\n',encoding='utf-8')
    (HERE/'requirements-python-audit-delta.txt').write_text(
      '# Delta relative to work/repro_v3_environment/requirements-python.lock.txt.\n'+
      ('# No additional packages: the baseline already contains rdata==1.1.0.\n' if not added else '')+
      '\n'.join(current[k]['line'] for k in added)+'\n',encoding='utf-8')
    (HERE/'requirements-python-optional-audit.txt').write_text(
      '# Explicit optional cross-language RDS review pin; already in the baseline.\n'+
      '# Install only on top of the matching baseline, without updating its dependencies.\n'+
      current['rdata']['line']+'\n',encoding='utf-8')
    tree=ast.parse(SCRIPT.read_text(encoding='utf-8'));imports=set()
    for node in ast.walk(tree):
        if isinstance(node,ast.Import):imports.update(x.name.split('.')[0] for x in node.names)
        elif isinstance(node,ast.ImportFrom) and node.module:imports.add(node.module.split('.')[0])
    third=sorted(imports-set(sys.stdlib_module_names))
    assert set(third)=={'numpy','pandas','scipy','liana','rdata'},third
    modules={}
    for name in third:
        module=importlib.import_module(name);dist=md.distribution(name);key=canonicalize_name(name)
        assert current[key]['version']==dist.version
        imported=getattr(module,'__version__',None)
        if imported is not None:assert str(imported)==dist.version,(name,imported,dist.version)
        modules[name]={'distribution_version':dist.version,'imported_version':imported,
          'module_path':str(Path(module.__file__).resolve()),'module_sha256':digest(module.__file__)}
    queue=[(name,()) for name in third];visited=set();records={};edges=[]
    while queue:
        name,extras=queue.pop();name=canonicalize_name(name);state=(name,tuple(sorted(extras)))
        if state in visited:continue
        visited.add(state);dist=md.distribution(name)
        assert name in current and current[name]['version']==dist.version,name
        metadata_text=dist.read_text('METADATA')
        records[name]={'name':dist.metadata['Name'],'version':dist.version,
          'requires_python':dist.metadata.get('Requires-Python'),'requires_dist':dist.requires or [],
          'metadata_sha256':hashlib.sha256(metadata_text.encode('utf-8')).hexdigest(),
          'home_page':dist.metadata.get('Home-page'),'project_urls':dist.metadata.get_all('Project-URL') or [],
          'metadata_path':str(getattr(dist,'_path',''))}
        for text in dist.requires or []:
            req=Requirement(text)
            active=req.marker is None or any(req.marker.evaluate({'extra':x}) for x in ('',)+tuple(extras))
            if not active:continue
            dep=canonicalize_name(req.name);actual=md.version(dep)
            assert dep in current and current[dep]['version']==actual,(name,dep)
            assert not req.specifier or req.specifier.contains(actual,prereleases=True),(name,dep,text,actual)
            edges.append({'from':name,'to':dep,'requirement':text,'installed_version':actual,'extras':sorted(req.extras)})
            queue.append((dep,tuple(sorted(req.extras))))
    check=subprocess.run([sys.executable,'-m','pip','check'],capture_output=True,text=True,encoding='utf-8',env=env)
    assert check.returncode==0,(check.stdout,check.stderr)
    rdata=md.distribution('rdata')
    (HERE/'rdata-1.1.0.METADATA.txt').write_text(rdata.read_text('METADATA'),encoding='utf-8')
    write(HERE/'python_dependency_closure.json',{'script':SCRIPT.relative_to(ROOT).as_posix(),
      'script_sha256':digest(SCRIPT),'third_party_imports':third,'packages':records,'active_dependency_edges':edges,
      'closure_packages':len(records),'closure_fully_pinned':True,
      'marker_scope':'Actual Python 3.12.14/Windows environment; propagated requested package extras; development/documentation extras omitted unless requested by a production dependency.'})
    report={'status':'current_environment_and_optional_audit_dependencies_verified',
      'captured_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),
      'python_executable':sys.executable,'platform':platform.platform(),'sys_path':sys.path,
      'baseline_lock_path':BASE.relative_to(ROOT).as_posix(),'baseline_lock_sha256':digest(BASE),
      'baseline_pins':len(baseline),'current_pins':len(current),'added':{k:current[k]['version'] for k in added},
      'rdata_already_present_in_baseline':baseline.get('rdata',{}).get('version')==current['rdata']['version'],
      'removed':removed,'changed':changed,'excluded':excluded,'script_imports':modules,
      'independent_review_required_dependency_closure_packages':len(records),
      'all_active_dependency_constraints_and_current_pins_satisfied':True,
      'pip_check_exit_code':check.returncode,'pip_check_stdout':check.stdout,'pip_check_stderr':check.stderr,
      'fresh_python_environment_restore_tested':False,'conda_solve_attempted':False,
      'baseline_lock_modified':False,'package_installed_or_updated':False,'model_fit_performed':False,
      'prior_full_independent_review_passed':'work/repro_v3_cellchat_independent_review.json',
      'prior_review_audit_sha256':digest(SCRIPT.with_suffix('.json'))}
    write(HERE/'python_current_audit.json',report)
    print(json.dumps({'status':report['status'],'baseline_pins':len(baseline),'current_pins':len(current),
      'added':report['added'],'changed':changed,'dependency_closure_packages':len(records)},indent=2))

if __name__=='__main__':main()
