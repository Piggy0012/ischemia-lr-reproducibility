"""Metadata-only CellChat resource overlap and dense-allocation feasibility audit."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'repro_v3_cellchat_probe'
def table(name):
    return pd.read_csv(OUT/f'CellChatDB.mouse_{name}.tsv',sep='\t',keep_default_na=False)
lr=table('interaction')
cx=table('complex').set_index('resource_row_id')
cf=table('cofactor').set_index('resource_row_id')
def expand(x,tab):
    if not x:return []
    if x in tab.index:
        return [str(v) for v in tab.loc[x].to_numpy() if str(v)]
    return [x]
def canonical(xs):return '_'.join(sorted(set(xs)))
protein=lr[lr.annotation!='Non-protein Signaling'].copy()
protein['canonical_ligand']=protein.ligand.map(lambda x:canonical(expand(x,cx)))
protein['canonical_receptor']=protein.receptor.map(lambda x:canonical(expand(x,cx)))
protein['n_rows_per_canonical_pair']=protein.groupby(['canonical_ligand','canonical_receptor']).interaction_name.transform('size')
old=pd.read_csv(ROOT/'literature/mouseconsensus.csv')
oldset={(canonical(str(a).split('_')),canonical(str(b).split('_'))) for a,b in zip(old.source_genesymbol,old.target_genesymbol)}
protein['in_liana_mouseconsensus']=[(a,b) in oldset for a,b in zip(protein.canonical_ligand,protein.canonical_receptor)]
protein['has_cofactor']=protein[['agonist','antagonist','co_A_receptor','co_I_receptor']].ne('').any(axis=1)
protein.to_csv(OUT/'CellChatDB.mouse_protein_canonical_overlap.tsv',sep='\t',index=False)
def genes_for(row,include_cofactors=True):
    genes=set(expand(row.ligand,cx)+expand(row.receptor,cx))
    if include_cofactors:
        for f in ['agonist','antagonist','co_A_receptor','co_I_receptor']:genes.update(expand(row[f],cf))
    return genes
protein_genes=set().union(*(genes_for(r) for _,r in protein.iterrows()))
rows=[]
covered_sets=[]
for acc in ['GSE174574','GSE245386']:
    for ap in sorted((ROOT/'repro_liana_rank_diagnostic/raw'/acc).glob('primary__GSM*.json')):
        a=json.loads(ap.read_text());sample=a['sample']
        genes=set(pd.read_csv(ROOT/f'revision_cache/{acc}/{sample}/genes.tsv.gz',sep='\t').symbol.astype(str))
        n_signal=len(protein_genes&genes);cells=a['n_context_cells'];k=len(a['context_counts'])
        covered=set(protein.loc[[genes_for(r).issubset(genes) for _,r in protein.iterrows()],'interaction_name'])
        covered_sets.append(covered)
        rows.append({'dataset':acc,'sample':sample,'condition':a['condition'],'context_cells':cells,'context_types':k,'native_protein_signaling_genes_in_matrix':n_signal,
                     'covered_native_protein_rows_all_LR_and_cofactor_genes':len(covered),'one_double_dense_signal_matrix_bytes':8*cells*n_signal,
                     'four_dense_copies_plus_boot_array_bytes':32*cells*n_signal+8*100*n_signal*k,
                     'two_probability_arrays_bytes':16*k*k*len(protein),
                     'estimate_note':'allocation accounting only; excludes R package baseline, original sparse input and temporary per-gene vectors; not measured peak RSS'})
cost=pd.DataFrame(rows);cost.to_csv(OUT/'sample_memory_estimates.tsv',sep='\t',index=False)
common=set.intersection(*covered_sets)
protein['covered_all_11_LR_cofactor_genes']=protein.interaction_name.isin(common)
protein.to_csv(OUT/'CellChatDB.mouse_protein_canonical_overlap.tsv',sep='\t',index=False)
common_unique=protein[(protein.n_rows_per_canonical_pair==1)&protein.in_liana_mouseconsensus&protein.covered_all_11_LR_cofactor_genes]
summary={'native_all_rows':len(lr),'native_protein_rows':len(protein),'native_unique_canonical_protein_pairs':len(protein[['canonical_ligand','canonical_receptor']].drop_duplicates()),
         'native_protein_rows_with_cofactor':int(protein.has_cofactor.sum()),'native_protein_genes_including_cofactors':len(protein_genes),
         'liana_mouseconsensus_rows':len(old),'liana_mouseconsensus_unique_canonical_pairs':len(oldset),
         'native_protein_rows_with_pair_in_mouseconsensus':int(protein.in_liana_mouseconsensus.sum()),
         'native_unique_protein_pairs_in_mouseconsensus':len(protein.loc[protein.in_liana_mouseconsensus,['canonical_ligand','canonical_receptor']].drop_duplicates()),
         'canonical_pairs_with_multiple_native_rows':len(protein.loc[protein.n_rows_per_canonical_pair>1,['canonical_ligand','canonical_receptor']].drop_duplicates()),
         'all_11_gene_covered_native_protein_rows':len(common),
         'all_11_gene_covered_unambiguous_shared_rows':len(common_unique),
         'all_11_gene_covered_unambiguous_shared_no_cofactor_rows':int((~common_unique.has_cofactor).sum()),
         'largest_allocation_estimate':cost.loc[cost.four_dense_copies_plus_boot_array_bytes.idxmax()].to_dict(),
         'original_samples':len(cost),'no_expression_matrices_loaded':True,'scope':'Gene coverage only, not detection eligibility or observed concordance.'}
(OUT/'resource_memory_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
