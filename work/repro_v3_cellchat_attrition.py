"""Explain every original primary target candidate's retention or loss; metadata only."""
from pathlib import Path
import hashlib,json
import pandas as pd
ROOT=Path(__file__).resolve().parent;BASE=ROOT/'repro_v3_cellchat';OUT=BASE/'tables';PROBE=ROOT/'repro_v3_cellchat_probe'
def read(p,**kw):return pd.read_csv(p,sep='\t',keep_default_na=False,**kw)
def canon(x):return '_'.join(sorted(set(str(x).split('_'))))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
oldpath=ROOT.parent/'outputs/reproducibility_v2/tables/null_fixed_candidate_keys.tsv'
old=read(oldpath);old=old[old.config=='primary'].rename(columns={'sender':'source','receiver':'target'}).copy()
assert len(old)==187
native=read(PROBE/'CellChatDB.mouse_interaction.tsv')
cx=read(PROBE/'CellChatDB.mouse_complex.tsv').set_index('resource_row_id');cf=read(PROBE/'CellChatDB.mouse_cofactor.tsv').set_index('resource_row_id')
def parts(x,tab):return [v for v in tab.loc[x].to_numpy() if v] if x in tab.index else ([x] if x else [])
native['canonical_ligand']=native.ligand.map(lambda x:canon('_'.join(parts(x,cx))))
native['canonical_receptor']=native.receptor.map(lambda x:canon('_'.join(parts(x,cx))))
protein=native[native.annotation!='Non-protein Signaling']
official_symbols=set(read(PROBE/'CellChatDB.mouse_geneInfo.tsv').Symbol)
protein_parts=[(r.resource_row_id,set(parts(r.ligand,cx)),set(parts(r.receptor,cx))) for r in protein.itertuples(index=False)]
control=read(BASE/'controlled_resource.tsv')
controlset=set(zip(control.canonical_ligand,control.canonical_receptor))
keys=read(OUT/'fixed_complete_target_keys.tsv');keyset=set(keys.itertuples(index=False,name=None))
native_pairs={(a,b) for a,b in zip(native.canonical_ligand,native.canonical_receptor)}
foldpairs={(a.lower(),b.lower()) for a,b in native_pairs}
core_genes={}
for p in sorted((ROOT/'revision_cache').glob('GSE*/GSM*/genes.tsv.gz')):
    if p.parts[-3] in ['GSE174574','GSE245386']:
        core_genes[p.parts[-2]]=set(read(p).symbol)
rows=[]
for r in old.itertuples(index=False):
    lig,rec=canon(r.ligand),canon(r.receptor);pair=(lig,rec)
    allmatches=native[(native.canonical_ligand==lig)&(native.canonical_receptor==rec)]
    pm=allmatches[allmatches.annotation!='Non-protein Signaling']
    kept=pair in controlset
    failures={}
    if len(pm)==1:
        rr=pm.iloc[0];needs=set(parts(rr.ligand,cx)+parts(rr.receptor,cx))
        for col in ['agonist','antagonist','co_A_receptor','co_I_receptor']:needs.update(parts(rr[col],cf))
        failures={s:sorted(needs-g) for s,g in core_genes.items() if needs-g}
    targetkept=(r.source,r.target,lig,rec) in keyset
    subunit_supersets=[name for name,ll,rr in protein_parts if set(lig.split('_'))<=ll and set(rec.split('_'))<=rr and (set(lig.split('_'))!=ll or set(rec.split('_'))!=rr)]
    reason=('retained' if targetkept else 'no_exact_complete_subunit_pair_in_native_CellChatDB' if not len(allmatches) else
            'native_match_only_nonprotein' if not len(pm) else 'multiple_native_rows_for_same_canonical_pair' if len(pm)>1 else
            'missing_LR_or_native_cofactor_gene_in_one_or_more_libraries' if failures else
            'not_in_frozen_controlled_resource' if not kept else 'not_in_all_11_LIANA_full_network_intersection')
    rows.append({'source':r.source,'target':r.target,'original_ligand':r.ligand,'original_receptor':r.receptor,
      'canonical_ligand':lig,'canonical_receptor':rec,'complex_order_changed':lig!=r.ligand or rec!=r.receptor,
      'n_native_rows_exact_canonical':len(allmatches),'n_native_protein_rows_exact_canonical':len(pm),
      'native_interaction_ids':'|'.join(pm.resource_row_id),'casefold_only_native_match':not len(allmatches) and (lig.lower(),rec.lower()) in foldpairs,
      'all_original_subunit_symbols_in_native_geneInfo':set(lig.split('_')+rec.split('_'))<=official_symbols,
      'native_rows_requiring_additional_subunits':'|'.join(subunit_supersets),
      'missing_genes_by_sample':json.dumps(failures,sort_keys=True),'in_controlled_1548':kept,'in_final_32':targetkept,'reason':reason})
d=pd.DataFrame(rows)
assert int(d.in_final_32.sum())==32 and set(zip(d.loc[d.in_final_32,'source'],d.loc[d.in_final_32,'target'],d.loc[d.in_final_32,'canonical_ligand'],d.loc[d.in_final_32,'canonical_receptor']))==keyset
assert not d.casefold_only_native_match.any(), 'Investigate case-only mismatches before interpreting resource loss'
assert not (d.in_controlled_1548 & ~d.in_final_32).any(),'Original complete targets should not vanish at full-network intersection'
d.to_csv(OUT/'original_187_to_controlled_32_attrition.tsv',sep='\t',index=False)
summary={'status':'complete','original_primary_target_candidates':187,'native_all_categories_exact_canonical_matches':int((d.n_native_rows_exact_canonical>0).sum()),
 'native_protein_exact_canonical_matches':int((d.n_native_protein_rows_exact_canonical>0).sum()),'native_protein_unique_row_matches':int((d.n_native_protein_rows_exact_canonical==1).sum()),
 'controlled_shared_1548_resource_matches':int(d.in_controlled_1548.sum()),'after_all_11_full_network_intersection':int(d.in_final_32.sum()),
 'additional_target_loss_at_full_network_intersection':int((d.in_controlled_1548 & ~d.in_final_32).sum()),'casefold_only_unmatched_pairs':int(d.casefold_only_native_match.sum()),
 'canonical_complex_order_changes':int(d.complex_order_changed.sum()),'reasons':d.reason.value_counts().to_dict(),
 'original_candidates_with_all_genes_as_exact_native_official_symbols':int(d.all_original_subunit_symbols_in_native_geneInfo.sum()),
 'excluded_candidates_with_native_superset_complex_match':int((~d.in_final_32 & d.native_rows_requiring_additional_subunits.ne('')).sum()),
 'resource_loss_interpretation':'Exact ligand/receptor subunit definitions differ. Do not collapse a native obligate receptor complex to a single component merely to force matches; no alias conversion or orthology mapping was applied.',
 'old_target_source_sha256':sha(oldpath),'controlled_resource_sha256':sha(BASE/'controlled_resource.tsv'),'native_resource_snapshot_sha256':sha(PROBE/'data/CellChatDB.mouse.rda'),
 'scope':'Resource and coverage audit, not an observed CellChat reproducibility result.'}
(OUT/'original_187_to_controlled_32_attrition.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
