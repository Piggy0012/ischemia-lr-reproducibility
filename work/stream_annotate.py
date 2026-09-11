"""Memory-bounded raw UMI QC and transparent lineage-marker annotation.
This script does not test disease associations or choose markers from DE results.
"""
from pathlib import Path
import gzip,json,sys,gc
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
MARKERS={
 'Astrocyte':['Aldh1l1','Slc1a3','Slc1a2','Gja1','Glul','Sox9','Aldoc'],
 'Endothelial':['Pecam1','Cdh5','Esam','Tek','Erg','Vwf','Kdr'],
 'Pericyte':['Pdgfrb','Rgs5','Kcnj8','Abcc9','Cspg4','Des'],
 'Smooth_muscle':['Acta2','Tagln','Myh11','Cnn1','Myl9'],
 'Fibroblast':['Col1a1','Col1a2','Dcn','Lum','Col3a1'],
 'Microglia':['P2ry12','Tmem119','Cx3cr1','Hexb','Csf1r'],
 'Myeloid':['Lyz2','Tyrobp','Aif1','Fcgr3','Ctss','C1qa'],
 'Oligodendrocyte':['Mbp','Plp1','Mog','Mag','Cnp','Mobp'],
 'OPC':['Pdgfra','Olig1','Olig2','Sox10','Gpr17'],
 'Neuron':['Snap25','Syt1','Rbfox3','Slc17a7','Gad1','Gad2'],
 'Ependymal':['Foxj1','Tppp3','Pifo','Dynlrb2','Ccdc153'],
 'Choroid':['Ttr','Krt18','Krt8','Aqp1','Folr1'],
 'Lymphocyte':['Cd3d','Cd3e','Trbc2','Cd79a','Ms4a1','Nkg7'],
 'Neutrophil':['S100a8','S100a9','Retnlg','Ly6g','Csf3r'],
 'Erythroid':['Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Alas2'],
}
def chunks(path):
    with gzip.open(path,'rt') as f:
        line=f.readline()
        while line.startswith('%'):line=f.readline()
        dims=list(map(int,line.split()))
        yield dims
        for chunk in pd.read_csv(f,sep=r'\s+',header=None,names=['g','c','x'],dtype=np.int32,chunksize=250000):
            a=chunk.to_numpy();a[:,:2]-=1;yield a
def process(acc):
    raw=ROOT/'data'/acc/'raw';out=ROOT/'processed'/acc;out.mkdir(parents=True,exist_ok=True)
    (out/'marker_panel.json').write_text(json.dumps(MARKERS,indent=2),encoding='utf-8')
    summary=[]
    for matrix in sorted(raw.glob('*matrix.mtx.gz')):
        prefix=matrix.name.removesuffix('_matrix.mtx.gz');sample=prefix.split('_')[0]
        dst=out/(sample+'_cells.tsv.gz')
        if dst.exists():
            q=pd.read_csv(dst,sep='\t'); print('EXISTS',sample,flush=True)
        else:
            print('ANNOTATE',sample,flush=True)
            feature=next(iter(list(raw.glob(prefix+'_genes.tsv.gz'))+list(raw.glob(prefix+'_features.tsv.gz'))))
            features=pd.read_csv(feature,sep='\t',header=None);symbols=features.iloc[:,1].astype(str).to_numpy()
            barcodes=pd.read_csv(raw/(prefix+'_barcodes.tsv.gz'),sep='\t',header=None).iloc[:,0].to_numpy()
            it=chunks(matrix);ng,nc,nnz=next(it)
            assert ng==len(symbols) and nc==len(barcodes)
            mg=sorted(set(sum(MARKERS.values(),[])));idx={g:i for i,g in enumerate(mg)}
            gl=np.array([idx.get(s,-1) for s in symbols]);mt=np.array([s.lower().startswith('mt-') for s in symbols])
            total=np.zeros(nc,np.float64);nfeat=np.zeros(nc,np.int32);mito=np.zeros(nc,np.float64);m=np.zeros((nc,len(mg)),np.float32)
            seen=0
            for a in it:
                assert a[:,2].min()>0
                g,c,x=a.T
                total+=np.bincount(c,weights=x,minlength=nc)
                nfeat+=np.bincount(c,minlength=nc).astype(np.int32)
                yes=mt[g];mito+=np.bincount(c[yes],weights=x[yes],minlength=nc)
                yes=gl[g]>=0;np.add.at(m,(c[yes],gl[g[yes]]),x[yes]);seen+=len(a)
            assert seen==nnz
            q=pd.DataFrame({'barcode':barcodes,'sample':sample,'prefix':prefix,'n_umis':total.astype(np.int64),'n_genes':nfeat,'pct_mt':mito/np.maximum(total,1)*100})
            q['qc_pass']=(q.n_genes.between(200,6000)&(q.n_umis>=500)&(q.pct_mt<=20))
            norm=np.log1p(m/np.maximum(total[:,None],1)*1e4)
            good=q.qc_pass.to_numpy();mu=norm[good].mean(axis=0);sd=norm[good].std(axis=0)
            z=np.clip((norm-mu)/np.maximum(sd,.25),-3,10)
            scores=[];detect=[]
            for typ,genes in MARKERS.items():
                ids=[idx[g] for g in genes];s=z[:,ids].mean(axis=1);d=(m[:,ids]>0).sum(axis=1)
                q['score_'+typ]=s;q['nmarkers_'+typ]=d;scores.append(s);detect.append(d)
            scores=np.array(scores).T;detect=np.array(detect).T;order=np.argsort(scores,axis=1);best=order[:,-1];second=order[:,-2]
            q['marker_type']=np.array(list(MARKERS))[best];q['identity_score']=scores[np.arange(nc),best]
            q['identity_margin']=q.identity_score.to_numpy()-scores[np.arange(nc),second]
            q['identity_nmarkers']=detect[np.arange(nc),best]
            q['identity_pass']=(q.identity_score>=.5)&(q.identity_margin>=.25)&(q.identity_nmarkers>=2)
            q['cell_type']=np.where(q.qc_pass&q.identity_pass,q.marker_type,'Unassigned')
            q['strict_identity']=(q.identity_score>=.75)&(q.identity_margin>=.5)&(q.identity_nmarkers>=3)
            q.to_csv(dst,sep='\t',index=False)
            features.iloc[:,:2].set_axis(['gene_id','symbol'],axis=1).to_csv(out/(sample+'_features.tsv.gz'),sep='\t',index=False)
            del m,norm,z,scores,detect;gc.collect()
        counts=q.cell_type.value_counts().to_dict()
        rec={'sample':sample,'n_input':len(q),'n_qc':int(q.qc_pass.sum()),'n_assigned':int((q.cell_type!='Unassigned').sum()),**counts};summary.append(rec)
        print(json.dumps(rec),flush=True)
    pd.DataFrame(summary).fillna(0).to_csv(out/'qc_summary.tsv',sep='\t',index=False)
if __name__=='__main__':process(sys.argv[1] if len(sys.argv)>1 else 'GSE174574')
