from pathlib import Path
import requests,json,hashlib
ROOT=Path(__file__).resolve().parent
p=ROOT/'spatial_geometry';p.mkdir(exist_ok=True)
u='https://api.github.com/repos/pachterlab/SpatialFeatureExperiment/git/trees/main?recursive=1'
r=requests.get(u,timeout=60);r.raise_for_status();tree=r.json();(p/'repository_tree.json').write_text(json.dumps(tree),encoding='utf-8')
for item in tree['tree']:
    if 'visium_row_col' in item['path']:
        print(item['path'])
        if item['path'].endswith(('.rda','.rds','.RData','.csv')):
            url='https://raw.githubusercontent.com/pachterlab/SpatialFeatureExperiment/'+tree['sha']+'/'+item['path'];a=requests.get(url,timeout=60);a.raise_for_status();dest=p/Path(item['path']).name;dest.write_bytes(a.content)
            (p/'manifest.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(a.content).hexdigest(),'commit':tree['sha']},indent=2),encoding='utf-8')
