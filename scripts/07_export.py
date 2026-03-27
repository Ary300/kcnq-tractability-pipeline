from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import EXPORT_DIR


def main() -> None:
    top10 = pd.read_csv(EXPORT_DIR / "top10_final.csv")
    clinical = pd.read_csv(EXPORT_DIR / "clinical_variants_scored.csv")
    docking_top10 = pd.read_csv(EXPORT_DIR / "top10_docking.csv")
    docking_matrix = pd.read_csv(EXPORT_DIR / "docking_matrix.csv")
    stats_summary = pd.read_csv(EXPORT_DIR / "section32_stats.csv")
    regression = pd.read_csv(EXPORT_DIR / "section32_logit.csv")
    sensitivity_grid = pd.read_csv(EXPORT_DIR / "weight_sensitivity_grid.csv") if (EXPORT_DIR / "weight_sensitivity_grid.csv").exists() else pd.DataFrame()
    sensitivity_candidates = pd.read_csv(EXPORT_DIR / "weight_sensitivity_candidates.csv") if (EXPORT_DIR / "weight_sensitivity_candidates.csv").exists() else pd.DataFrame()
    sensitivity_effects = pd.read_csv(EXPORT_DIR / "weight_sensitivity_effects.csv") if (EXPORT_DIR / "weight_sensitivity_effects.csv").exists() else pd.DataFrame()
    lit = pd.read_csv(EXPORT_DIR / "kcnq4_literature_comparison.csv")
    kcnq1_validation = pd.read_csv(EXPORT_DIR / "kcnq1_candidate_validation.csv") if (EXPORT_DIR / "kcnq1_candidate_validation.csv").exists() else pd.DataFrame()
    brewer_overlap = pd.read_csv(EXPORT_DIR / "kcnq1_brewer_overlap.csv") if (EXPORT_DIR / "kcnq1_brewer_overlap.csv").exists() else pd.DataFrame()
    brewer_mechanism = pd.read_csv(EXPORT_DIR / "kcnq1_brewer_mechanism_summary.csv") if (EXPORT_DIR / "kcnq1_brewer_mechanism_summary.csv").exists() else pd.DataFrame()
    brewer_mechanism_omnibus = pd.read_csv(EXPORT_DIR / "kcnq1_brewer_mechanism_omnibus.csv") if (EXPORT_DIR / "kcnq1_brewer_mechanism_omnibus.csv").exists() else pd.DataFrame()
    ml277_status = pd.read_csv(EXPORT_DIR / "kcnq1_ml277_status.csv") if (EXPORT_DIR / "kcnq1_ml277_status.csv").exists() else pd.DataFrame()
    conservation = pd.read_csv(EXPORT_DIR / "top10_conservation.csv") if (EXPORT_DIR / "top10_conservation.csv").exists() else pd.DataFrame()
    am_concordance_summary = pd.read_csv(EXPORT_DIR / "alphamissense_concordance_summary.csv") if (EXPORT_DIR / "alphamissense_concordance_summary.csv").exists() else pd.DataFrame()
    am_auc = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_alphamissense_auc.csv") if (EXPORT_DIR / "kcnq2_bhatt_alphamissense_auc.csv").exists() else pd.DataFrame()
    am_high_subset_auc = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_alphamissense_high_subset_auc.csv") if (EXPORT_DIR / "kcnq2_bhatt_alphamissense_high_subset_auc.csv").exists() else pd.DataFrame()
    revel_auc = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_revel_auc.csv") if (EXPORT_DIR / "kcnq2_bhatt_revel_auc.csv").exists() else pd.DataFrame()
    bhatt_overlap = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_overlap.csv") if (EXPORT_DIR / "kcnq2_bhatt_overlap.csv").exists() else pd.DataFrame()
    bhatt_summary = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_summary.csv") if (EXPORT_DIR / "kcnq2_bhatt_summary.csv").exists() else pd.DataFrame()
    bhatt_retig_overlap = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_retigabine_overlap.csv") if (EXPORT_DIR / "kcnq2_bhatt_retigabine_overlap.csv").exists() else pd.DataFrame()
    bhatt_retig_summary = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_retigabine_summary.csv") if (EXPORT_DIR / "kcnq2_bhatt_retigabine_summary.csv").exists() else pd.DataFrame()
    bhatt_mechanism_summary = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_mechanism_summary.csv") if (EXPORT_DIR / "kcnq2_bhatt_mechanism_summary.csv").exists() else pd.DataFrame()
    bhatt_hetero_overlap = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_heteromeric_overlap.csv") if (EXPORT_DIR / "kcnq2_bhatt_heteromeric_overlap.csv").exists() else pd.DataFrame()
    bhatt_hetero_summary = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_heteromeric_summary.csv") if (EXPORT_DIR / "kcnq2_bhatt_heteromeric_summary.csv").exists() else pd.DataFrame()
    bhatt_hetero_mech = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_heteromeric_mechanism_summary.csv") if (EXPORT_DIR / "kcnq2_bhatt_heteromeric_mechanism_summary.csv").exists() else pd.DataFrame()
    bhatt_operating_metrics = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_operating_point_metrics.csv") if (EXPORT_DIR / "kcnq2_bhatt_operating_point_metrics.csv").exists() else pd.DataFrame()
    bhatt_confusion = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_confusion_matrix.csv") if (EXPORT_DIR / "kcnq2_bhatt_confusion_matrix.csv").exists() else pd.DataFrame()
    bhatt_delong = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_delong.csv") if (EXPORT_DIR / "kcnq2_bhatt_delong.csv").exists() else pd.DataFrame()
    bhatt_loo_auc = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_leave_one_out_auc.csv") if (EXPORT_DIR / "kcnq2_bhatt_leave_one_out_auc.csv").exists() else pd.DataFrame()
    vanoye_overlap = pd.read_csv(EXPORT_DIR / "kcnq1_vanoye_overlap.csv") if (EXPORT_DIR / "kcnq1_vanoye_overlap.csv").exists() else pd.DataFrame()
    vanoye_summary = pd.read_csv(EXPORT_DIR / "kcnq1_vanoye_summary.csv") if (EXPORT_DIR / "kcnq1_vanoye_summary.csv").exists() else pd.DataFrame()

    if not conservation.empty:
        top10 = top10.merge(conservation, on=["gene", "protein_change", "residue_num"], how="left")

    final_top10 = top10[
        [
            "final_rank" if "final_rank" in top10.columns else ("priority_rank" if "priority_rank" in top10.columns else "rescue_rank"),
            "gene",
            "protein_change",
            "residue_num",
            "region_label",
            "clinical_significance",
            "trait",
            "pocket_dist_A",
            "pocket_proximity_score",
            "path_score",
            "region_drug_score",
            "rarity_score",
            "structural_opportunity_score" if "structural_opportunity_score" in top10.columns else "rescue_score",
            "tractability_modifier",
            "tractability_class",
            "tractability_basis",
            "rescue_priority_score",
            "dG_WT",
            "docking_wt_score",
            "final_rescue_score",
            "rescue_tier",
            "candidate_drug",
            "drug_phase",
            "paralog_identity_label" if "paralog_identity_label" in top10.columns else "region_label",
            "cross_species_identity_label" if "cross_species_identity_label" in top10.columns else "rarity_score",
            "paralog_conservation_status" if "paralog_conservation_status" in top10.columns else "region_label",
            "cross_species_conservation_status" if "cross_species_conservation_status" in top10.columns else "region_label",
            "consurf_like_grade_1to9" if "consurf_like_grade_1to9" in top10.columns else "rarity_score",
            "structural_coords_available",
            "source",
        ]
    ].copy()
    final_top10.columns = [
        "Rank",
        "Gene",
        "Variant",
        "Residue",
        "Structural Region",
        "ClinVar Significance",
        "Disease/Trait",
        "Pocket Distance (Å)",
        "Pocket Proximity Score",
        "Pathogenicity Score",
        "Region Drug Score",
        "Rarity Score",
        "Structural Opportunity Score",
        "Tractability Modifier",
        "Tractability Class",
        "Tractability Basis",
        "Rescue Priority Score",
        "WT Docking ΔG (kcal/mol)",
        "WT Docking Score",
        "Final Rescue Score",
        "Rescue Tier",
        "Candidate Drug",
        "Drug Phase",
        "Paralog Identity (X/5)",
        "Cross-Species Identity (X/4)",
        "Paralog Conservation Status",
        "Cross-Species Conservation Status",
        "ConSurf-Like Grade (1-9)",
        "Structural Coords Available",
        "Source",
    ]

    clinical_out = clinical[
        [
            "gene",
            "protein_change",
            "residue_num",
            "clinical_significance",
            "trait",
            "region_label",
            "pocket_dist_A",
            "pocket_proximity_score",
            "structural_coords_available",
            "clinvar_tier_score",
            "cadd_phred",
            "cadd_norm",
            "alphamissense_score",
            "path_score_imputed",
            "path_score",
            "region_drug_score",
            "rarity_score",
            "rescue_score",
            "structural_opportunity_score",
            "tractability_modifier",
            "tractability_class",
            "tractability_basis",
            "rescue_priority_score",
            "gnomad_af",
            "source",
            "hgvs_c",
            "variant_id",
        ]
    ].sort_values("rescue_score", ascending=False)

    final_top10.to_csv(EXPORT_DIR / "kcnq_top10_rescue_candidates_final.csv", index=False)
    clinical_out.to_csv(EXPORT_DIR / "kcnq_all_annotated_variants_final.csv", index=False)

    excel_path = EXPORT_DIR / "kcnq_rescue_analysis_final.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        final_top10.to_excel(writer, sheet_name="Top10_Rescue_Candidates", index=False)
        clinical_out.to_excel(writer, sheet_name="All_Annotated_Variants", index=False)
        docking_top10.to_excel(writer, sheet_name="Docking_Top10", index=False)
        docking_matrix.to_excel(writer, sheet_name="Docking_Matrix", index=False)
        stats_summary.to_excel(writer, sheet_name="Section3_2_Stats", index=False)
        regression.to_excel(writer, sheet_name="Section3_2_Logit", index=False)
        sensitivity_grid.to_excel(writer, sheet_name="Sensitivity_Grid", index=False)
        sensitivity_candidates.to_excel(writer, sheet_name="Sensitivity_Candidates", index=False)
        sensitivity_effects.to_excel(writer, sheet_name="Sensitivity_Effects", index=False)
        lit.to_excel(writer, sheet_name="KCNQ4_Literature", index=False)
        kcnq1_validation.to_excel(writer, sheet_name="KCNQ1_Candidate_Valid", index=False)
        brewer_overlap.to_excel(writer, sheet_name="KCNQ1_Brewer_Overlap", index=False)
        brewer_mechanism.to_excel(writer, sheet_name="KCNQ1_Brewer_Mech", index=False)
        brewer_mechanism_omnibus.to_excel(writer, sheet_name="KCNQ1_Brewer_Omnibus", index=False)
        ml277_status.to_excel(writer, sheet_name="KCNQ1_ML277_Status", index=False)
        conservation.to_excel(writer, sheet_name="Top10_Conservation", index=False)
        am_concordance_summary.to_excel(writer, sheet_name="AM_Concordance_Sum", index=False)
        am_auc.to_excel(writer, sheet_name="KCNQ2_AM_vs_Struct", index=False)
        am_high_subset_auc.to_excel(writer, sheet_name="KCNQ2_AM_HighSubset", index=False)
        revel_auc.to_excel(writer, sheet_name="KCNQ2_REVEL_AUC", index=False)
        bhatt_overlap.to_excel(writer, sheet_name="KCNQ2_Bhatt_Overlap", index=False)
        bhatt_summary.to_excel(writer, sheet_name="KCNQ2_Bhatt_Summary", index=False)
        bhatt_retig_overlap.to_excel(writer, sheet_name="KCNQ2_Bhatt_Retig", index=False)
        bhatt_retig_summary.to_excel(writer, sheet_name="KCNQ2_Bhatt_Retig_Sum", index=False)
        bhatt_mechanism_summary.to_excel(writer, sheet_name="KCNQ2_Bhatt_Mech", index=False)
        bhatt_hetero_overlap.to_excel(writer, sheet_name="KCNQ2_Bhatt_Hetero", index=False)
        bhatt_hetero_summary.to_excel(writer, sheet_name="KCNQ2_Bhatt_HeteroSum", index=False)
        bhatt_hetero_mech.to_excel(writer, sheet_name="KCNQ2_Bhatt_HetMech", index=False)
        bhatt_operating_metrics.to_excel(writer, sheet_name="KCNQ2_Bhatt_OpPoint", index=False)
        bhatt_confusion.to_excel(writer, sheet_name="KCNQ2_Bhatt_ConfMat", index=False)
        bhatt_delong.to_excel(writer, sheet_name="KCNQ2_Bhatt_DeLong", index=False)
        bhatt_loo_auc.to_excel(writer, sheet_name="KCNQ2_Bhatt_LOO_AUC", index=False)
        vanoye_overlap.to_excel(writer, sheet_name="KCNQ1_Vanoye_Overlap", index=False)
        vanoye_summary.to_excel(writer, sheet_name="KCNQ1_Vanoye_Summary", index=False)

    print(excel_path)


if __name__ == "__main__":
    main()
