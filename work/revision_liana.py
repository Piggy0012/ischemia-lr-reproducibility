"""Actual five-method LIANA ensemble, separately for each animal library."""
from pathlib import Path
import os,sys,json,gc,hashlib,importlib.metadata
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
import numpy as np,pandas as pd,anndata as ad,liana as li
from revision_data import load_sample,samples
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'revision_results/liana'

def process(acc,sample,config):
    folder=OUT/acc;folder.mkdir(parents=True,exist_ok=True)
    dest=folder/f'{config}__{sample}.tsv.gz';audit=folder/f'{config}__{sample}.json'
    if os.environ.get('ISCHEMIA_REVISION_FORCE')!='1' and dest.exists() and audit.exists() and json.loads(audit.read_text())['parameters']['return_all_lrs'] is False:
        print('LIANA EXISTS',acc,sample,config,flush=True);return
    x,q,genes=load_sample(acc,sample);keep=(q.cell_type!='Unassigned').to_numpy()
    if config=='reference_singlet':
        r=pd.read_csv(ROOT/f'revision_results/identity/{acc}/{sample}_identity.tsv.gz',sep='\t')
        assert np.array_equal(q.barcode,r.barcode)
        target=q.cell_type.isin(['Astrocyte','Endothelial']).to_numpy()
        supported=(r.whole_brain_broad==q.cell_type).to_numpy()&(r.whole_brain_probability.to_numpy()>=.5)
        keep&=(~r.predicted_doublet.to_numpy(bool))&((~target)|supported)
    else:assert config=='primary'
    groupcounts=q.loc[keep,'cell_type'].value_counts()
    retained=groupcounts[groupcounts>=30].index
    keep&=q.cell_type.isin(retained).to_numpy()
    assert all(t in retained for t in ['Astrocyte','Endothelial'])
    sub=q.loc[keep];y=x[keep].astype(np.float32);del x;gc.collect()
    y=y.multiply((1e4/sub.n_umis.to_numpy()).astype(np.float32)[:,None]).tocsr();y.data=np.log1p(y.data)
    a=ad.AnnData(y,obs=pd.DataFrame({'cell_type':pd.Categorical(sub.cell_type.to_numpy())},index=sub.barcode.astype(str)),var=pd.DataFrame(index=genes))
    resource_path=ROOT/'literature/mouseconsensus.csv'
    resource=pd.read_csv(resource_path)[['source_genesymbol','target_genesymbol']].drop_duplicates().rename(columns={'source_genesymbol':'ligand','target_genesymbol':'receptor'})
    print('LIANA START',acc,sample,config,'shape',a.shape,'groups',sub.cell_type.value_counts().to_dict(),flush=True)
    # All cell-type pairs are computed: restricting pairs upstream also alters
    # expression context and NATMI denominators. Only exports are restricted.
    result=li.mt.rank_aggregate(a,groupby='cell_type',resource=resource,expr_prop=.1,min_cells=30,
        use_raw=False,n_perms=100,seed=20260910,n_jobs=1,return_all_lrs=False,inplace=False,verbose=False)
    mask=((result.source=='Astrocyte')&(result.target=='Endothelial'))|((result.source=='Endothelial')&(result.target=='Astrocyte'))
    target_result=result.loc[mask].copy()
    target_result.insert(0,'sample',sample);target_result.insert(0,'dataset',acc);target_result.insert(0,'config',config)
    condition='Sham' if ('sham' in q.prefix.iloc[0].lower() or '_WTC' in q.prefix.iloc[0]) else 'MCAO'
    target_result['condition']=condition
    target_result['expression_eligible']=True
    assert not target_result.duplicated(['source','target','ligand_complex','receptor_complex']).any()
    target_result.to_csv(dest,sep='\t',index=False)
    report={'dataset':acc,'sample':sample,'condition':condition,'config':config,
      'liana_version':importlib.metadata.version('liana'),'methods':[m.method_name for m in li.mt.rank_aggregate.methods],
      'parameters':{'expr_prop':.1,'min_cells':30,'n_perms':100,'seed':20260910,'n_jobs':1,'return_all_lrs':False,'use_raw':False},
      'n_input_cells':len(q),'n_context_cells':a.n_obs,'n_context_genes':a.n_vars,'context_counts':sub.cell_type.value_counts().to_dict(),
      'computed_all_cell_type_pairs':True,'all_pair_rows':len(result),'exported_target_rows':len(target_result),
      'resource_sha256':hashlib.sha256(resource_path.read_bytes()).hexdigest(),'columns':list(target_result.columns),
      'interpretation':'CellPhoneDB permutation P values are within-library cell-label statistics, not animal-level disease P values. Consensus ranks are relative priorities in the full computed network, not measured signaling probabilities. Default return_all_lrs=False reports expression-eligible interactions only. Missing rows are not observed zero scores and must not be filled with zero for animal-level contrasts.'}
    audit.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('LIANA DONE',sample,config,'all',len(result),'target',len(target_result),flush=True)
    del result,target_result,a,y;gc.collect()

if __name__=='__main__':
    # No arguments runs all 22 library/config combinations sequentially.
    accs=[sys.argv[1]] if len(sys.argv)>1 else ['GSE174574','GSE245386']
    for acc in accs:
        for sample in ([sys.argv[2]] if len(sys.argv)>2 else samples(acc)):
            for config in ([sys.argv[3]] if len(sys.argv)>3 else ['primary','reference_singlet']):process(acc,sample,config)
