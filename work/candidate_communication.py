from pathlib import Path
import sys,json,itertools
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parent
def process(acc,config='primary'):
    folder=ROOT/'processed'/acc;out=ROOT/'results'/acc;out.mkdir(parents=True,exist_ok=True)
    manifest=pd.read_csv(folder/'sample_cell_counts.tsv',sep='\t');manifest=manifest[manifest.config==config]
    samples=sorted(manifest['sample'].unique());condition=manifest.drop_duplicates('sample').set_index('sample').loc[samples,'condition'];case=(condition=='MCAO').to_numpy()
    means={};fractions={};de={}
    for typ in ['Astrocyte','Endothelial']:
        means[typ]=pd.concat([pd.read_csv(folder/(s+'_mean_logcp10k.tsv.gz'),sep='\t',index_col=0)[config+'__'+typ].rename(s) for s in samples],axis=1).fillna(0)
        fractions[typ]=pd.concat([pd.read_csv(folder/(s+'_fractions.tsv.gz'),sep='\t',index_col=0)[config+'__'+typ].rename(s) for s in samples],axis=1).fillna(0)
        de[typ]=pd.read_csv(out/f'{config}__{typ}__de.tsv.gz',sep='\t',index_col=0)
    lr=pd.read_csv(ROOT/'literature/mouseconsensus.csv').drop_duplicates(['source_genesymbol','target_genesymbol'])
    records=[];score_records=[]
    for sender,receiver in [('Astrocyte','Endothelial'),('Endothelial','Astrocyte')]:
        for row in lr.itertuples(index=False):
            ligand=row.source_genesymbol;receptor=row.target_genesymbol;lg=ligand.split('_');rg=receptor.split('_')
            if not all(g in means[sender].index for g in lg) or not all(g in means[receiver].index for g in rg):continue
            ld=fractions[sender].loc[lg].min(axis=0).to_numpy();rd=fractions[receiver].loc[rg].min(axis=0).to_numpy()
            eligible=[]
            for threshold in [.05,.1,.2]:
                lpass=max((ld[case]>=threshold).sum(),(ld[~case]>=threshold).sum())>=2
                rpass=max((rd[case]>=threshold).sum(),(rd[~case]>=threshold).sum())>=2
                eligible.append(bool(lpass and rpass))
            if not eligible[0]:continue
            l=means[sender].loc[lg].min(axis=0).to_numpy();r=means[receiver].loc[rg].min(axis=0).to_numpy()
            score=np.sqrt(l*r);delta=float(score[case].mean()-score[~case].mean())
            ttest=stats.ttest_ind(score[case],score[~case],equal_var=False)
            perm=[]
            for chosen in itertools.combinations(range(len(samples)),int(case.sum())):
                mask=np.zeros(len(samples),bool);mask[list(chosen)]=True;perm.append(score[mask].mean()-score[~mask].mean())
            exact=float((np.abs(perm)>=abs(delta)-1e-12).mean())
            lf=de[sender].reindex(lg);rf=de[receiver].reindex(rg)
            l_changed=bool(((lf.padj<.05)&(lf.log2FoldChange.abs()>=.5)).any());r_changed=bool(((rf.padj<.05)&(rf.log2FoldChange.abs()>=.5)).any())
            rec={'dataset':acc,'config':config,'sender':sender,'receiver':receiver,'ligand':ligand,'receptor':receptor,'eligible_05':eligible[0],'eligible_10':eligible[1],'eligible_20':eligible[2],'mean_control_score':float(score[~case].mean()),'mean_stroke_score':float(score[case].mean()),'score_difference':delta,'welch_p':float(ttest.pvalue),'exact_permutation_p':exact,'ligand_mean_log2fc':float(lf.log2FoldChange.mean()),'receptor_mean_log2fc':float(rf.log2FoldChange.mean()),'ligand_min_padj':float(lf.padj.min()),'receptor_min_padj':float(rf.padj.min()),'ligand_any_changed':l_changed,'receptor_any_changed':r_changed,'candidate_changed_component':l_changed or r_changed,'ligand_min_detection_stroke':float(ld[case].mean()),'receptor_min_detection_stroke':float(rd[case].mean())};records.append(rec)
            for i,s in enumerate(samples):score_records.append({'sender':sender,'receiver':receiver,'ligand':ligand,'receptor':receptor,'sample':s,'condition':condition.loc[s],'score':float(score[i]),'ligand_logexpr':float(l[i]),'receptor_logexpr':float(r[i])})
    table=pd.DataFrame(records);table['welch_fdr']=np.nan
    for gate in ['eligible_05','eligible_10','eligible_20']:
        ix=table[gate];table.loc[ix,'welch_fdr_'+gate]=multipletests(table.loc[ix,'welch_p'].fillna(1),method='fdr_bh')[1]
    table.to_csv(out/f'{config}__communication.tsv.gz',sep='\t',index=False)
    pd.DataFrame(score_records).to_csv(out/f'{config}__communication_sample_scores.tsv.gz',sep='\t',index=False)
    primary=table[table.eligible_10];print(acc,config,'eligible',len(primary),'changed component',int(primary.candidate_changed_component.sum()),'WelchFDR05',int((primary.welch_fdr_eligible_10<.05).sum()),'minExactP',primary.exact_permutation_p.min(),flush=True)
if __name__=='__main__':process(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'primary')
