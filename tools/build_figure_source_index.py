"""Explicit figure provenance mapping; index existing bytes without plotting."""
from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parent.parent
OUT = 'outputs/reproducibility_v2/'
T = OUT + 'tables/'
S = OUT + 'figure_source_data/'


def main():
    definitions = [
        ('Figure1','Figure1_samples_and_marker_groups','work/repro_core_figures.py',
         ['work/results/identity_marker_audit.tsv'] + [f'work/processed/{a}/sample_cell_counts.tsv' for a in ['GSE174574','GSE245386']],
         [S+'figure1_animal_counts.tsv',S+'figure1_markers_GSE174574.tsv',S+'figure1_markers_GSE245386.tsv']),
        ('Figure2','Figure2_expression_priority_and_null','work/repro_core_figures.py',
         [T+n for n in ['null_fixed_common_observed_effects.tsv.gz','null_fixed_common_summary.tsv','null_dynamic_gate_all_allocations.tsv.gz','null_dynamic_gate_summary.tsv']],
         [S+n for n in ['figure2_common_effects.tsv.gz','figure2_common_null_summary.tsv','figure2_original_gate_null.tsv']]),
        ('Figure3','Figure3_implementation_and_rank_universe','work/repro_rank_figure.py',
         [T+n for n in ['rank_background_summary.tsv','rank_background_context.tsv','rank_background_zero_effects.tsv']],
         [S+n for n in ['figure3_comparison_and_denominators.tsv','figure3_global_universe_context.tsv','figure3_zero_effect_sources.tsv']]),
        ('Figure4','Figure4_ambient_sensitivity','work/repro_ambient_figure.py',
         [T+n for n in ['ambient_markers_per_animal.tsv','ambient_rank_matched_summary.tsv','ambient_concordance_summary.tsv','null_dynamic_gate_summary.tsv']],
         [S+n for n in ['figure4_paired_markers.tsv','figure4_matched_rank_comparisons.tsv','figure4_gate_coverage.tsv']]),
        ('Figure5','Figure5_third_cohort_transfer','work/repro_third_figure.py',
         ['work/repro_third_cohort/library_qc_summary.tsv',T+'third_rank_diagnostic_summary.tsv'],
         [S+n for n in ['figure5_library_counts.tsv','figure5_fixed_17_library_comparison.tsv','figure5_all_metric_and_universe_results.tsv']]),
        ('FigureS1','FigureS1_reference_identity','work/revision_identity_figures.py',[],
         ['outputs/revised/figure_source_data/'+n for n in ['figure6_sample_metadata.tsv','figure6_target_proportions.tsv','figure6_full_label_confusion.tsv','figure6_coarse_confusion_by_sample.tsv','figure6_confusion_animal_mean.tsv']]),
        ('FigureS2','FigureS2_doublet_distributions','work/revision_identity_figures.py',[],
         ['outputs/revised/figure_source_data/figureS3_histogram_bins.tsv','outputs/revised/figure_source_data/figureS3_sample_doublet_audits.tsv']),
        ('FigureS3','FigureS3_barrier_selection','work/revision_barrier_figure.py',
         ['outputs/revised/tables/barrier_genes_cell_selection_sensitivity.tsv','work/gene_panels.json'],
         ['outputs/revised/figure_source_data/figureS4_barrier_effects.tsv']),
        ('FigureS4','FigureS4_sorted_effects_and_intervals','work/repro_sorted_descriptive.py',
         ['work/sorted_rna/astro_endothelial_normalized.tsv.gz','work/sorted_rna/diff_expression.tsv.gz'],
         [S+'sorted_fixed_components_per_pool.tsv',T+'sorted_fixed_component_effects_ci.tsv']),
        ('FigureS5','FigureS5_spatial_expression_context','work/make_figures.py',
         [f'work/results/spatial/{sample}_spots.tsv.gz' for sample in ['GSM7437221','GSM7437222']],
         [f'outputs/figure_source_data/figure5_{sample}.tsv.gz' for sample in ['GSM7437221','GSM7437222']]),
    ]
    rows = []
    for fig, base, script, inputs, sources in definitions:
        copy_script = 'work/repro_copy_supplement_figures.py' if fig in ['FigureS1','FigureS2','FigureS3','FigureS5'] else ''
        def add(role,path,note=''):
            p = ROOT / path
            assert p.is_file(), path
            rows.append({'figure_id':fig,'figure_basename':base,'role':role,'path':path,
                         'plot_script':script,'copy_script':copy_script,'bytes':p.stat().st_size,
                         'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'note':note})
        add('plot_script',script)
        if copy_script:
            add('versioned_export_copy_script',copy_script)
            add('versioned_export_provenance',OUT+'supplement_figure_provenance.json',
                'Previously verified vector/PNG export copied unchanged into v2 numbering; source location is retained in this audit.')
        for path in inputs:
            add('actual_plot_input',path)
        for path in sources:
            add('exported_figure_source_data',path)
        for ext in ['png','pdf','svg']:
            add('figure_export',OUT+'figures/'+base+'.'+ext)
        if fig in ['FigureS1','FigureS2']:
            add('actual_plot_input','work/revision_sample_manifest.tsv')
            for acc in ['GSE174574','GSE245386']:
                add('actual_plot_input',f'work/processed/{acc}/qc_summary.tsv')
                for suffix in ['_audit.json','_identity.tsv.gz','_simulated_doublet_scores.tsv.gz']:
                    samples = range(5319987,5319993) if acc=='GSE174574' else range(7841720,7841725)
                    paths = [ROOT / 'work/revision_results/identity' / acc / (f'GSM{s}'+suffix) for s in samples]
                    for p in paths:
                        add('actual_plot_input',p.relative_to(ROOT).as_posix(),
                            'Shared complete-input loader reads labels, real/simulated doublet scores and audits; no new model fit.')
        if fig == 'Figure1':
            for acc in ['GSE174574','GSE245386']:
                add('upstream_cell_identity_definition',f'work/processed/{acc}/marker_panel.json',
                    'Upstream frozen assignment panel; the displayed Figure 1 marker subset is defined in the plotting script.')
        if fig == 'Figure5':
            add('upstream_cell_identity_definition','work/repro_third_cohort/processed/GSE332910/marker_panel.json')
            for name in ['third_rank_diagnostic_original_pair_on99.tsv','third_rank_diagnostic_original_pair_on99_effects.tsv.gz','third_rank_diagnostic_original_pair_on99_audit.json']:
                add('companion_analysis_not_plotted',T+name,'Original two cohorts restricted to the same 99-candidate set; distinct from original-versus-third contrasts plotted in Figure 5.')
        if fig == 'FigureS4':
            add('complete_unselected_descriptive_results',T+'sorted_all_gene_descriptive_effects.tsv.gz',
                '49,914 gene/cell-type paired-pool effects with pointwise intervals; fixed illustrative genes only are displayed.')
            add('descriptive_source_audit',T+'sorted_descriptive_audit.json')
    target = ROOT / 'FIGURE_SOURCE_INDEX.tsv'
    s6root='outputs/reproducibility_v3_cellchat/'
    if (ROOT/s6root/'figures/FigureS6_cellchat_controlled_extension.png').exists():
        fig='FigureS6';base='FigureS6_cellchat_controlled_extension'
        script='work/repro_v3_cellchat_figure.py';copy_script=''
        add('plot_script',script)
        for name in ['target_disease_effects.tsv.gz','cross_cohort_concordance.tsv','zero_strength_target_audit.tsv','tie_zero_audit.tsv','summary_audit.json']:
            add('actual_plot_input','work/repro_v3_cellchat/tables/'+name)
        for name in ['figureS6_candidate_effects.tsv','figureS6_concordance.tsv','figureS6_zero_tie_floors.tsv']:
            add('exported_figure_source_data',s6root+'figure_source_data/'+name)
        for ext in ['png','pdf','svg']:
            add('figure_export',s6root+'figures/'+base+'.'+ext)
        add('visual_review',s6root+'figure_source_data/figureS6_visual_review.json')
    with target.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    print(json.dumps({'figures':len({r['figure_id'] for r in rows}),'indexed_source_and_export_rows':len(rows),'all_paths_exist':True},indent=2))


if __name__=='__main__':
    main()
