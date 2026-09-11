from pathlib import Path
import pickle,json,pandas as pd
p=Path('work/identity_reference_review');out=[]
for f in p.glob('*.pkl'):
 with f.open('rb') as h:m=pickle.load(h)
 print(f.name,type(m),list(m),flush=True)
 clf=m['Model'];scaler=m['Scaler_'];genes=set(clf.features)
 d={'model':f.name,'classes':list(clf.classes_),'n_features':len(genes),'gene_head':list(clf.features)[:15],'coef_shape':list(clf.coef_.shape),'scaler_with_mean':scaler.with_mean,'description':m['description']}
 matches=[]
 for acc in ['GSE174574','GSE245386']:
  for q in sorted((Path('work/data')/acc/'raw').glob('*genes.tsv.gz'))+sorted((Path('work/data')/acc/'raw').glob('*features.tsv.gz')):
   g=pd.read_csv(q,sep='\t',header=None).iloc[:,1].astype(str)
   common=len(set(g)&genes);matches.append({'file':q.name,'n_query_unique':len(set(g)),'model_genes_matched':common,'model_gene_fraction':common/len(genes)})
 d['gene_matches']=matches;out.append(d);print(json.dumps(d,default=str),flush=True)
(p/'model_inspection.json').write_text(json.dumps(out,indent=2,default=str))

