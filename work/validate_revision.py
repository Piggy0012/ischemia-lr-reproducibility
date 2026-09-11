"""Numerical and provenance checks; not additional biological validation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/revised'
checks=[]
def check(ok,name):
    if not ok:raise AssertionError(name)
    checks.append(name)
summary=json.loads((OUT/'tables/summary.json').read_text())
check(summary['status']=='complete','All revision inputs complete')
check(not summary['missing_scores_imputed'],'LIANA scores not zero imputed')
check(not summary['liana_internal_pvalues_used_for_disease_inference'],'Animal inference distinct from internal cell-label P values')
manifest=pd.read_csv(ROOT/'revision_sample_manifest.tsv',sep='\t')
for acc in ['GSE174574','GSE245386']:
    count=pd.read_csv(ROOT/f'revision_analysis/processed/{acc}/sample_cell_counts.tsv',sep='\t')
    for row in manifest[manifest.dataset==acc].itertuples(index=False):
        sample=row.sample;stem=ROOT/f'revision_results/identity/{acc}/{sample}'
        r=pd.read_csv(str(stem)+'_identity.tsv.gz',sep='\t');a=json.loads(Path(str(stem)+'_audit.json').read_text())
        check(r.barcode.is_unique and len(r)==a['n_qc'],sample+' unique cells and QC count')
        check(int(r.predicted_doublet.sum())==a['scrublet']['n_predicted'],sample+' doublet count')
        check(np.array_equal(r.doublet_score>a['scrublet']['threshold'],r.predicted_doublet),sample+' saved threshold reproduces flags')
        for model in a['models'].values():check(model['compatibility']['labels_identical'] and model['compatibility']['probabilities_identical'],sample+' model '+model['model']+' compatibility')
        singlet=~r.predicted_doublet
        for config in ['singlet','reference_singlet','reference_only']:
            for typ in ['Astrocyte','Endothelial']:
                original=r.cell_type==typ;ref=(r.whole_brain_broad==typ)&(r.whole_brain_probability>=.5)
                mask=singlet&(original if config=='singlet' else (original&ref if config=='reference_singlet' else ref))
                c=count[(count['sample']==sample)&(count.config==config)&(count.cell_type==typ)]
                check(len(c)==1 and mask.sum()==c.n_cells.iloc[0] and mask.sum()>=30,sample+' '+config+' '+typ+' cells')
                check(int(r.loc[mask,'n_umis'].sum())==int(c.total_target_counts.iloc[0]),sample+' '+config+' '+typ+' total all-gene counts')
        for config in ['primary','reference_singlet']:
            la=json.loads((ROOT/f'revision_results/liana/{acc}/{config}__{sample}.json').read_text())
            check(la['computed_all_cell_type_pairs'] and la['parameters']['return_all_lrs'] is False,sample+' '+config+' LIANA final settings')
effects=pd.read_csv(OUT/'tables/liana_animal_level_effects.tsv.gz',sep='\t')
check(effects.animal_exact_permutation_p.ge(.1-1e-12).all(),'All exact animal P values respect design resolution')
check(effects.animal_welch_fdr.dropna().between(0,1).all(),'Finite animal BH values bounded')
examples=pd.read_csv(OUT/'tables/liana_fixed_examples_animal_effects.tsv',sep='\t')
check(examples.loc[examples.status=='not_complete_case','score_difference'].isna().all(),'Incomplete fixed examples have no imputed disease effect')
for dataset,n in [('GSE174574',6),('GSE245386',5)]:
    check(effects[effects.dataset==dataset].n_sham.eq(3).all() and effects[effects.dataset==dataset].n_mcao.eq(n-3).all(),dataset+' biological replicate counts')
for folder in ['tables','figure_source_data']:
    check(any((OUT/folder).iterdir()),folder+' present')
for f in ROOT.glob('revision_*.py'):compile(f.read_text(encoding='utf-8'),str(f),'exec')
report={'status':'PASS','checks_passed':len(checks),'checks':checks,'inference_unit':'animal library',
        'scope':'Checks cell counts, summed counts, stored doublet decisions, model compatibility, final LIANA settings, missingness and animal inference. It does not establish biological identity or functional communication.',
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'validation_revision.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('Revision validation PASS',len(checks))
