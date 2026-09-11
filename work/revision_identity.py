"""Independent mouse-reference annotation and per-library Scrublet audit."""
from pathlib import Path
import os,sys,json,gc,logging
ROOT=Path(__file__).resolve().parent
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
os.environ['CELLTYPIST_FOLDER']=str(ROOT/'revision_cache/celltypist')
import numpy as np,pandas as pd
from scipy.special import expit
import anndata as ad,scanpy as sc,celltypist
from revision_data import load_sample,samples

OUT=ROOT/'revision_results/identity';OUT.mkdir(parents=True,exist_ok=True)
celltypist.logger.set_level(logging.WARNING)
MODELS={
 'whole_brain':ROOT/'identity_reference_review/Mouse_Whole_Brain.pkl',
 'cortex_hippocampus':ROOT/'identity_reference_review/Mouse_Isocortex_Hippocampus.pkl'}

def broad(label,reference):
    if reference=='cortex_hippocampus':
        return {'Astro':'Astrocyte','Endo':'Endothelial','Micro-PVM':'Microglia/Macrophage','Oligo':'Oligodendrocyte','SMC-Peri':'Mural','VLMC':'Fibroblast'}.get(label,'Neuron/Other')
    code=int(label.split(' ')[0])
    if code in [317,318,319,320]:return 'Astrocyte'
    return {321:'Astroependymal',322:'Tanycyte',323:'Ependymal',324:'Hypendymal',325:'Choroid',326:'OPC',327:'Oligodendrocyte',328:'OEC',329:'ABC',330:'Fibroblast',331:'Pericyte',332:'Smooth_muscle',333:'Endothelial',334:'Microglia',335:'Macrophage',336:'Monocytes',337:'Dendritic',338:'Lymphocyte'}.get(code,'Neuron/Other')

def annotate(x,q,genes,model_path):
    model=celltypist.models.Model.load(model_path.as_posix())
    matched=np.flatnonzero(genes.isin(model.features));mi=pd.Index(model.features).get_indexer(genes[matched])
    labels=[];probs=[];compat={}
    for start in range(0,len(q),256):
        stop=min(start+256,len(q));y=x[start:stop].astype(np.float64)
        y=y.multiply((1e4/q.n_umis.to_numpy()[start:stop])[:,None]).tocsr();y.data=np.log1p(y.data)
        a=ad.AnnData(y,obs=pd.DataFrame(index=q.barcode.iloc[start:stop].astype(str)),var=pd.DataFrame(index=genes))
        res=celltypist.annotate(a,model=model,mode='best match',majority_voting=False)
        label=res.predicted_labels.predicted_labels.astype(str).to_numpy();prob=res.probability_matrix.max(axis=1).to_numpy()
        if start==0:
            z=(y[:,matched].toarray()-(model.scaler.mean_[mi] if model.scaler.with_mean else 0))/model.scaler.scale_[mi]
            z[z>10]=10;manual=z@model.classifier.coef_[:,mi].T+model.classifier.intercept_
            assert np.allclose(manual,res.decision_matrix.to_numpy(),rtol=1e-8,atol=1e-8)
            assert np.array_equal(model.classifier.classes_[manual.argmax(axis=1)],label)
            assert np.allclose(expit(manual).max(axis=1),prob,rtol=1e-8,atol=1e-8)
            compat={'checked_cells':len(label),'decision_max_abs_error':float(np.max(abs(manual-res.decision_matrix.to_numpy()))),'labels_identical':True,'probabilities_identical':True}
        labels.extend(label);probs.extend(prob)
        del res,a,y
    return np.asarray(labels),np.asarray(probs),{'model':model_path.name,'matched_features':len(mi),'model_features':len(model.features),'compatibility':compat,'majority_voting':False,'probability_note':'One-versus-rest sigmoid score, not a calibrated posterior or ground truth'}

def process(acc,sample):
    folder=OUT/acc;folder.mkdir(exist_ok=True);dest=folder/f'{sample}_identity.tsv.gz';auditfile=folder/f'{sample}_audit.json'
    force=os.environ.get('ISCHEMIA_REVISION_FORCE')=='1'
    if auditfile.exists() and not force:print('IDENTITY EXISTS',acc,sample,flush=True);return
    x,q,genes=load_sample(acc,sample);audits={'dataset':acc,'sample':sample,'n_qc':len(q),'models':{}}
    if dest.exists() and not force:
        r=pd.read_csv(dest,sep='\t');assert np.array_equal(q.barcode,r.barcode)
        audits=json.loads((folder/f'{sample}_models_audit.json').read_text())
    else:
        r=q[['barcode','sample','cell_type','n_umis','n_genes','pct_mt']].copy()
        for name,path in MODELS.items():
            print('REFERENCE',acc,sample,name,flush=True)
            labels,prob,report=annotate(x,q,genes,path)
            r[name+'_label']=labels;r[name+'_broad']=[broad(v,name) for v in labels];r[name+'_probability']=prob;audits['models'][name]=report
        r.to_csv(dest,sep='\t',index=False)
        (folder/f'{sample}_models_audit.json').write_text(json.dumps(audits,indent=2),encoding='utf-8')
    print('SCRUBLET',acc,sample,flush=True)
    a=ad.AnnData(x,obs=pd.DataFrame(index=q.barcode.astype(str)),var=pd.DataFrame(index=genes))
    sc.pp.scrublet(a,expected_doublet_rate=.05,sim_doublet_ratio=2.,n_prin_comps=30,use_approx_neighbors=False,random_state=20260910,verbose=True)
    r['doublet_score']=a.obs.doublet_score.to_numpy();r['predicted_doublet']=a.obs.predicted_doublet.to_numpy()
    s=a.uns['scrublet'];threshold=float(s['threshold']);sim=np.asarray(s['doublet_scores_sim'])
    pd.DataFrame({'simulated_score':sim}).to_csv(folder/f'{sample}_simulated_doublet_scores.tsv.gz',sep='\t',index=False)
    r.to_csv(dest,sep='\t',index=False)
    audits['scrublet']={'expected_doublet_rate':.05,'sim_doublet_ratio':2.,'n_prin_comps':30,'random_state':20260910,'threshold':threshold,'n_simulated':len(sim),'n_predicted':int(r.predicted_doublet.sum()),'predicted_fraction':float(r.predicted_doublet.mean()),'automatic_threshold':True,'ambient_RNA_corrected':False}
    audits['rule_target_counts']=r.cell_type.value_counts().to_dict()
    audits['reference_target_counts']=r.whole_brain_broad.value_counts().to_dict()
    auditfile.write_text(json.dumps(audits,indent=2),encoding='utf-8')
    print('IDENTITY DONE',sample,'doublets',audits['scrublet']['n_predicted'],'threshold',threshold,flush=True)
    del a,x,r;gc.collect()

if __name__=='__main__':
    accs=[sys.argv[1]] if len(sys.argv)>1 else ['GSE174574','GSE245386']
    for acc in accs:
        for sample in ([sys.argv[2]] if len(sys.argv)>2 else samples(acc)):process(acc,sample)
