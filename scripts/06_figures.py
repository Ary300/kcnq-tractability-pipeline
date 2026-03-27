from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import EXPORT_DIR
from kcnq_pipeline.figures import (
    save_alphamissense_concordance_figure,
    save_benchmark_figure,
    save_docking_bar_chart,
    save_docking_heatmap,
    save_modular_benchmark_figure,
    save_operating_point_confusion_figure,
    save_pose_render_figure,
    save_sensitivity_heatmap,
    save_structural_annotation_figure,
    save_top10_figure,
    save_umap,
)


def main() -> None:
    clinical = pd.read_csv(EXPORT_DIR / "clinical_variants_scored.csv")
    top10 = pd.read_csv(EXPORT_DIR / "top10_scored.csv")
    docking = pd.read_csv(EXPORT_DIR / "docking_matrix.csv")
    docking_top10 = pd.read_csv(EXPORT_DIR / "top10_docking.csv")
    stats_summary = pd.read_csv(EXPORT_DIR / "section32_stats.csv") if (EXPORT_DIR / "section32_stats.csv").exists() else pd.DataFrame()
    merged_top10 = top10.merge(
        docking_top10[["gene", "protein_change", "dG_WT", "docking_wt_score"]],
        on=["gene", "protein_change"],
        how="left",
    )
    merged_top10["final_rescue_score"] = 0.80 * merged_top10["rescue_priority_score"] + 0.20 * merged_top10["docking_wt_score"].fillna(0.0)
    merged_top10["rescue_tier"] = merged_top10["final_rescue_score"].map(
        lambda s: "HIGH" if s >= 0.85 else ("MODERATE" if s >= 0.75 else "LOW")
    )
    docking["variant_label"] = docking["gene"] + " " + docking["protein_change"]
    save_top10_figure(merged_top10)
    save_umap(clinical, merged_top10)
    save_structural_annotation_figure(clinical, stats_summary)
    save_docking_bar_chart(docking_top10)
    save_docking_heatmap(docking)
    kcnq2_overlap = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_retigabine_overlap.csv") if (EXPORT_DIR / "kcnq2_bhatt_retigabine_overlap.csv").exists() else pd.DataFrame()
    kcnq1_overlap = pd.read_csv(EXPORT_DIR / "kcnq1_vanoye_overlap.csv") if (EXPORT_DIR / "kcnq1_vanoye_overlap.csv").exists() else pd.DataFrame()
    sensitivity = pd.read_csv(EXPORT_DIR / "weight_sensitivity_effects.csv") if (EXPORT_DIR / "weight_sensitivity_effects.csv").exists() else pd.DataFrame()
    am_concordance = pd.read_csv(EXPORT_DIR / "alphamissense_concordance_scatter.csv") if (EXPORT_DIR / "alphamissense_concordance_scatter.csv").exists() else pd.DataFrame()
    am_summary = pd.read_csv(EXPORT_DIR / "alphamissense_concordance_summary.csv") if (EXPORT_DIR / "alphamissense_concordance_summary.csv").exists() else pd.DataFrame()
    am_auc = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_alphamissense_auc.csv") if (EXPORT_DIR / "kcnq2_bhatt_alphamissense_auc.csv").exists() else pd.DataFrame()
    modular_auc = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_modular_predictor_auc.csv") if (EXPORT_DIR / "kcnq2_bhatt_modular_predictor_auc.csv").exists() else pd.DataFrame()
    op_metrics = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_operating_point_metrics.csv") if (EXPORT_DIR / "kcnq2_bhatt_operating_point_metrics.csv").exists() else pd.DataFrame()
    op_confusion = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_confusion_matrix.csv") if (EXPORT_DIR / "kcnq2_bhatt_confusion_matrix.csv").exists() else pd.DataFrame()
    save_benchmark_figure(kcnq2_overlap, kcnq1_overlap, merged_top10)
    save_alphamissense_concordance_figure(am_concordance, am_summary, am_auc)
    save_modular_benchmark_figure(kcnq2_overlap, modular_auc)
    save_operating_point_confusion_figure(op_metrics, op_confusion)
    save_pose_render_figure()
    save_sensitivity_heatmap(sensitivity)
    merged_top10.to_csv(EXPORT_DIR / "top10_final.csv", index=False)
    print("Figures written")


if __name__ == "__main__":
    main()
