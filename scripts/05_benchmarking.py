from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.benchmarking import (
    benchmark_alphamissense_high_subsets,
    benchmark_alphamissense_vs_structural,
    benchmark_bhatt_delong,
    benchmark_bhatt_heteromeric_rescue,
    benchmark_bhatt_mechanism_strata,
    benchmark_bhatt_operating_point,
    benchmark_brewer_kcnq1_mechanisms,
    benchmark_bhatt_retigabine,
    benchmark_overlap,
    benchmark_revel_vs_structural,
    build_kcnq1_candidate_validation,
    build_kcnq4_literature_table,
    build_ml277_status_table,
    compute_alphamissense_concordance,
    compute_top10_paralog_conservation,
    leave_one_component_out_auc,
    parse_brewer_kcnq1_mave,
    parse_bhatt_kcnq2,
    parse_vanoye_kcnq1,
    section32_stats,
    weight_sensitivity,
)
from kcnq_pipeline.config import BENCHMARK_DIR, EXPORT_DIR


def main() -> None:
    clinical = pd.read_csv(EXPORT_DIR / "clinical_variants_scored.csv")
    top10 = pd.read_csv(EXPORT_DIR / "top10_scored.csv")
    stats_summary, regression = section32_stats(clinical)
    sensitivity_grid, sensitivity_candidates, sensitivity_effects = weight_sensitivity(clinical)
    conservation_path = EXPORT_DIR / "top10_conservation.csv"
    if conservation_path.exists():
        conservation = pd.read_csv(conservation_path)
    else:
        conservation = compute_top10_paralog_conservation(top10)
    am_concordance, am_concordance_summary = compute_alphamissense_concordance(clinical)
    bhatt = parse_bhatt_kcnq2() if (BENCHMARK_DIR / "bhatt_kcnq2_supp_table4.xlsx").exists() else pd.DataFrame()
    vanoye = parse_vanoye_kcnq1() if (BENCHMARK_DIR / "vanoye_kcnq1_dataset.xlsx").exists() else pd.DataFrame()
    brewer = parse_brewer_kcnq1_mave()
    kcnq4_table = build_kcnq4_literature_table(clinical, top10)
    kcnq1_validation = build_kcnq1_candidate_validation(top10, vanoye, brewer)
    brewer_overlap, brewer_mech_summary, brewer_mech_omnibus = benchmark_brewer_kcnq1_mechanisms(clinical, brewer)
    ml277_status = build_ml277_status_table()
    stats_summary.to_csv(EXPORT_DIR / "section32_stats.csv", index=False)
    regression.to_csv(EXPORT_DIR / "section32_logit.csv", index=False)
    sensitivity_grid.to_csv(EXPORT_DIR / "weight_sensitivity_grid.csv", index=False)
    sensitivity_candidates.to_csv(EXPORT_DIR / "weight_sensitivity_candidates.csv", index=False)
    sensitivity_effects.to_csv(EXPORT_DIR / "weight_sensitivity_effects.csv", index=False)
    kcnq4_table.to_csv(EXPORT_DIR / "kcnq4_literature_comparison.csv", index=False)
    kcnq1_validation.to_csv(EXPORT_DIR / "kcnq1_candidate_validation.csv", index=False)
    brewer_overlap.to_csv(EXPORT_DIR / "kcnq1_brewer_overlap.csv", index=False)
    brewer_mech_summary.to_csv(EXPORT_DIR / "kcnq1_brewer_mechanism_summary.csv", index=False)
    brewer_mech_omnibus.to_csv(EXPORT_DIR / "kcnq1_brewer_mechanism_omnibus.csv", index=False)
    ml277_status.to_csv(EXPORT_DIR / "kcnq1_ml277_status.csv", index=False)
    conservation.to_csv(EXPORT_DIR / "top10_conservation.csv", index=False)
    am_concordance.to_csv(EXPORT_DIR / "alphamissense_concordance_scatter.csv", index=False)
    am_concordance_summary.to_csv(EXPORT_DIR / "alphamissense_concordance_summary.csv", index=False)
    if not bhatt.empty:
        bhatt_overlap, bhatt_summary = benchmark_overlap(clinical, bhatt, gene="KCNQ2", benchmark_prefix="bhatt")
        bhatt_retig_overlap, bhatt_retig_summary = benchmark_bhatt_retigabine(clinical, bhatt)
        bhatt_mechanism_summary = benchmark_bhatt_mechanism_strata(bhatt_retig_overlap)
        bhatt_hetero_overlap, bhatt_hetero_summary, bhatt_hetero_mechanism = benchmark_bhatt_heteromeric_rescue(clinical, bhatt)
        bhatt_operating_metrics, bhatt_confusion = benchmark_bhatt_operating_point(clinical, bhatt)
        bhatt_delong = benchmark_bhatt_delong(clinical, bhatt)
        loo_auc = leave_one_component_out_auc(clinical, bhatt)
        am_auc = benchmark_alphamissense_vs_structural(clinical, bhatt)
        am_high_subset_auc = benchmark_alphamissense_high_subsets(clinical, bhatt)
        bhatt_overlap.to_csv(EXPORT_DIR / "kcnq2_bhatt_overlap.csv", index=False)
        bhatt_summary.to_csv(EXPORT_DIR / "kcnq2_bhatt_summary.csv", index=False)
        bhatt_retig_overlap.to_csv(EXPORT_DIR / "kcnq2_bhatt_retigabine_overlap.csv", index=False)
        bhatt_retig_summary.to_csv(EXPORT_DIR / "kcnq2_bhatt_retigabine_summary.csv", index=False)
        bhatt_mechanism_summary.to_csv(EXPORT_DIR / "kcnq2_bhatt_mechanism_summary.csv", index=False)
        bhatt_hetero_overlap.to_csv(EXPORT_DIR / "kcnq2_bhatt_heteromeric_overlap.csv", index=False)
        bhatt_hetero_summary.to_csv(EXPORT_DIR / "kcnq2_bhatt_heteromeric_summary.csv", index=False)
        bhatt_hetero_mechanism.to_csv(EXPORT_DIR / "kcnq2_bhatt_heteromeric_mechanism_summary.csv", index=False)
        bhatt_operating_metrics.to_csv(EXPORT_DIR / "kcnq2_bhatt_operating_point_metrics.csv", index=False)
        bhatt_confusion.to_csv(EXPORT_DIR / "kcnq2_bhatt_confusion_matrix.csv", index=False)
        bhatt_delong.to_csv(EXPORT_DIR / "kcnq2_bhatt_delong.csv", index=False)
        loo_auc.to_csv(EXPORT_DIR / "kcnq2_bhatt_leave_one_out_auc.csv", index=False)
        am_auc.to_csv(EXPORT_DIR / "kcnq2_bhatt_alphamissense_auc.csv", index=False)
        am_high_subset_auc.to_csv(EXPORT_DIR / "kcnq2_bhatt_alphamissense_high_subset_auc.csv", index=False)
        revel_scores_path = EXPORT_DIR / "kcnq2_bhatt_revel_scores.csv"
        if revel_scores_path.exists():
            revel_scores = pd.read_csv(revel_scores_path)
            revel_auc = benchmark_revel_vs_structural(clinical, bhatt, revel_scores)
            revel_auc.to_csv(EXPORT_DIR / "kcnq2_bhatt_revel_auc.csv", index=False)
    if not vanoye.empty:
        vanoye_overlap, vanoye_summary = benchmark_overlap(clinical, vanoye, gene="KCNQ1", benchmark_prefix="vanoye")
        vanoye_overlap.to_csv(EXPORT_DIR / "kcnq1_vanoye_overlap.csv", index=False)
        vanoye_summary.to_csv(EXPORT_DIR / "kcnq1_vanoye_summary.csv", index=False)
    print("Benchmark/stat tables written")


if __name__ == "__main__":
    main()
