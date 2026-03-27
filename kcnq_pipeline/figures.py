from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from Bio.PDB import PDBParser
from sklearn.metrics import auc, roc_curve
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats

from .config import FIGURES_DIR, ROOT


REGION_COLORS = {
    "Pore/Selectivity Filter": "#E63946",
    "Ligand-Binding Pocket": "#F4A261",
    "Pore Domain (S5-S6)": "#E76F51",
    "Voltage-Sensing Domain": "#457B9D",
    "CaM/Interface Region": "#2A9D8F",
    "Trafficking/Assembly": "#8338EC",
    "Unknown/Unresolved": "#ADB5BD",
}
TIER_COLORS = {"HIGH": "#E63946", "MODERATE": "#F4A261", "LOW": "#ADB5BD"}
CLINVAR_GROUP_COLORS = {
    "P/LP": "#E63946",
    "VUS/conflicting": "#457B9D",
    "Other": "#ADB5BD",
}
MECHANISM_COLORS = {
    "dominant_negative": "#8D99AE",
    "full_rescue": "#2A9D8F",
    "partial_rescue": "#E9C46A",
    "no_rescue_annotated": "#E76F51",
    "unassigned": "#ADB5BD",
}
MODEL_COLORS = {
    "alphamissense_only": "#1D4E89",
    "am_x_tractability": "#F28E2B",
    "rescue_priority_only": "#2A9D8F",
    "revel_only": "#111111",
    "revel_x_tractability": "#E63946",
    "structural_only": "#7B2CBF",
}


def save_top10_figure(df_top10: pd.DataFrame) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df = df_top10.copy().sort_values("final_rescue_score", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = [TIER_COLORS[t] for t in df["rescue_tier"]]
    labels = [f"{g} {v}" for g, v in zip(df["gene"], df["protein_change"])]
    ax.barh(labels, df["final_rescue_score"], color=colors)
    ax.set_xlabel("Final Rescue Score")
    ax.set_title("Top 10 Rescue Candidates")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig5_top10_ranked.png", dpi=150)
    plt.savefig(FIGURES_DIR / "fig5_top10_ranked.svg")
    plt.close(fig)


def save_umap(df_clinical: pd.DataFrame, df_top10: pd.DataFrame) -> None:
    use = df_clinical[df_clinical["residue_num"].notna()].copy()
    features = ["residue_num", "path_score", "pocket_proximity_score", "region_drug_score", "rarity_score"]
    X = StandardScaler().fit_transform(use[features].fillna(0).values)
    try:
        from umap import UMAP

        emb = UMAP(n_components=2, random_state=42, n_neighbors=20, min_dist=0.1).fit_transform(X)
        title = "Variant UMAP"
    except Exception:
        emb = PCA(n_components=2, random_state=42).fit_transform(X)
        title = "Variant Feature Projection (PCA fallback)"
    use["UMAP1"], use["UMAP2"] = emb[:, 0], emb[:, 1]
    top_keys = set(df_top10.apply(lambda r: f"{r['gene']} {r['protein_change']}", axis=1))
    use["is_top10"] = use.apply(lambda r: f"{r['gene']} {r['protein_change']}" in top_keys, axis=1)
    fig, ax = plt.subplots(figsize=(8, 6))
    for region, color in REGION_COLORS.items():
        mask = use["region_label"] == region
        if mask.any():
            ax.scatter(use.loc[mask, "UMAP1"], use.loc[mask, "UMAP2"], s=10, c=color, alpha=0.6, label=region)
    top = use[use["is_top10"]]
    ax.scatter(top["UMAP1"], top["UMAP2"], c="black", marker="*", s=90, label="Top 10")
    ax.legend(fontsize=7, loc="best")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig1_umap.png", dpi=150)
    plt.savefig(FIGURES_DIR / "fig1_umap.svg")
    plt.close(fig)


def save_docking_heatmap(matrix_df: pd.DataFrame) -> None:
    if matrix_df.empty:
        return
    pivot = matrix_df.pivot(index="variant_label", columns="drug_name", values="dG_WT")
    fig, ax = plt.subplots(figsize=(12, 7))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="viridis_r", ax=ax)
    ax.set_title("Docking Matrix (WT ΔG)")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig4_docking_matrix.png", dpi=150)
    plt.savefig(FIGURES_DIR / "fig4_docking_matrix.svg")
    plt.close(fig)


def save_structural_annotation_figure(df_clinical: pd.DataFrame, stats_summary: pd.DataFrame) -> None:
    tmp = df_clinical.copy()
    tmp["sig_group"] = np.where(
        tmp["clinical_significance"].isin(["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]),
        "P/LP",
        np.where(tmp["clinical_significance"] == "Uncertain significance", "VUS", "Other"),
    )
    use = tmp[tmp["sig_group"].isin(["P/LP", "VUS"])].copy()
    region_counts = (
        use.groupby(["sig_group", "region_label"]).size().reset_index(name="n")
        .pivot(index="sig_group", columns="region_label", values="n")
        .fillna(0)
    )
    region_counts = region_counts[[c for c in REGION_COLORS if c in region_counts.columns]]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    region_counts.plot(kind="bar", stacked=True, color=[REGION_COLORS[c] for c in region_counts.columns], ax=axes[0])
    axes[0].set_title("Structural Region Distribution")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Variant Count")
    axes[0].legend(fontsize=7, frameon=False)

    violin_df = use[use["pocket_dist_A"].notna() & use["sig_group"].isin(["P/LP", "VUS"])].copy()
    sns.violinplot(data=violin_df, x="sig_group", y="pocket_dist_A", palette={"P/LP": "#E76F51", "VUS": "#457B9D"}, ax=axes[1])
    axes[1].set_title("Pocket Distance by Clinical Class")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("Pocket Distance (Å)")
    if not stats_summary.empty:
        row = stats_summary.iloc[0]
        axes[1].text(
            0.03,
            0.97,
            f"p={row['p_value']:.3g}\nd={row['cohens_d']:.2f}\nΔmedian={row['median_difference']:.2f} Å",
            transform=axes[1].transAxes,
            va="top",
        )
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig2_structural_annotation.png", dpi=150)
    plt.savefig(FIGURES_DIR / "fig2_structural_annotation.svg")
    plt.close(fig)


def save_docking_bar_chart(top10_docking: pd.DataFrame) -> None:
    if top10_docking.empty:
        return
    df = top10_docking.copy().sort_values("dG_WT", ascending=True)
    df["variant_label"] = df["gene"] + " " + df["protein_change"]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(df["variant_label"], df["dG_WT"], color="#457B9D")
    ax.set_xlabel("WT Docking ΔG (kcal/mol)")
    ax.set_title("WT Docking Energies for Priority Candidates")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig3_wt_docking.png", dpi=150)
    plt.savefig(FIGURES_DIR / "fig3_wt_docking.svg")
    plt.close(fig)


def save_benchmark_figure(kcnq2_overlap: pd.DataFrame, kcnq1_overlap: pd.DataFrame, df_top10: pd.DataFrame | None = None) -> None:
    if kcnq2_overlap.empty or "bhatt_retig_current_rel_wt_ctrl" not in kcnq2_overlap.columns:
        return
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    scatter_ax = axes[0, 0]
    scatter_df = kcnq2_overlap.dropna(subset=["rescue_score", "bhatt_retig_current_rel_wt_ctrl"]).copy()
    for mechanism, color in MECHANISM_COLORS.items():
        sub = scatter_df[scatter_df["bhatt_mechanism_class"] == mechanism]
        if sub.empty:
            continue
        scatter_ax.scatter(
            sub["rescue_score"],
            sub["bhatt_retig_current_rel_wt_ctrl"],
            s=55,
            alpha=0.8,
            c=color,
            label=mechanism.replace("_", " "),
        )
    if len(scatter_df) >= 2:
        slope, intercept = np.polyfit(scatter_df["rescue_score"], scatter_df["bhatt_retig_current_rel_wt_ctrl"], 1)
        xs = np.linspace(scatter_df["rescue_score"].min(), scatter_df["rescue_score"].max(), 100)
        scatter_ax.plot(xs, slope * xs + intercept, color="black", linewidth=1)
        rho, p_value = stats.spearmanr(
            scatter_df["rescue_score"], scatter_df["bhatt_retig_current_rel_wt_ctrl"], nan_policy="omit"
        )
        scatter_ax.text(0.03, 0.96, f"structural rho={rho:.2f}\np={p_value:.3g}", transform=scatter_ax.transAxes, va="top")
    scatter_ax.set_title("KCNQ2 Retigabine Response by Mechanism")
    scatter_ax.set_xlabel("Structural Opportunity Score")
    scatter_ax.set_ylabel("Retigabine Current Relative to WT")
    scatter_ax.legend(frameon=False, fontsize=8, loc="lower left")

    roc_ax = axes[0, 1]
    roc_specs = [
        ("bhatt_retig_reaches_wt_current", "Strict WT-current endpoint", "rescue_score", "#E76F51", "-"),
        ("bhatt_retig_reaches_wt_current", "Strict WT-current endpoint", "rescue_priority_score", "#2A9D8F", "--"),
        ("bhatt_retig_responsive_any", "Responsive endpoint", "rescue_score", "#457B9D", "-"),
        ("bhatt_retig_responsive_any", "Responsive endpoint", "rescue_priority_score", "#264653", "--"),
    ]
    for endpoint_col, endpoint_label, score_col, color, linestyle in roc_specs:
        if score_col not in kcnq2_overlap.columns:
            continue
        roc_df = kcnq2_overlap.dropna(subset=[endpoint_col, score_col]).copy()
        if len(roc_df) < 3 or len(roc_df[endpoint_col].astype(int).unique()) < 2:
            continue
        fpr, tpr, _ = roc_curve(roc_df[endpoint_col].astype(int), roc_df[score_col])
        auc_val = auc(fpr, tpr)
        score_label = "priority" if score_col == "rescue_priority_score" else "structural"
        roc_ax.plot(fpr, tpr, color=color, linestyle=linestyle, linewidth=2, label=f"{endpoint_label} ({score_label}, AUC={auc_val:.2f})")
    roc_ax.plot([0, 1], [0, 1], color="#ADB5BD", linestyle=":", linewidth=1)
    roc_ax.set_title("ROC Under Two Endpoint Definitions")
    roc_ax.set_xlabel("False Positive Rate")
    roc_ax.set_ylabel("True Positive Rate")
    roc_ax.legend(frameon=False, fontsize=8, loc="lower right")

    box_ax = axes[1, 0]
    order = ["dominant_negative", "partial_rescue", "full_rescue", "no_rescue_annotated", "unassigned"]
    box_df = kcnq2_overlap.copy()
    box_df["mechanism_label"] = box_df["bhatt_mechanism_class"].str.replace("_", " ")
    sns.boxplot(
        data=box_df,
        x="bhatt_mechanism_class",
        y="rescue_score",
        order=order,
        palette=MECHANISM_COLORS,
        ax=box_ax,
    )
    box_ax.set_title("Structural Opportunity by Bhatt Mechanism Class")
    box_ax.set_xlabel("")
    box_ax.set_ylabel("Structural Opportunity Score")
    box_ax.set_xticklabels([label.replace("_", "\n") for label in order], rotation=0)

    table_ax = axes[1, 1]
    table_ax.axis("off")
    if df_top10 is None:
        df_top10 = pd.DataFrame()
    table_df = df_top10.copy()
    if not table_df.empty:
        table_df["Variant"] = table_df["gene"] + " " + table_df["protein_change"]
        use_cols = ["Variant", "candidate_drug", "tractability_class", "tractability_modifier", "rescue_priority_score"]
        use_cols = [col for col in use_cols if col in table_df.columns]
        table_df = table_df[use_cols].head(10).copy()
        rename_map = {
            "candidate_drug": "Drug",
            "tractability_class": "Tractability Class",
            "tractability_modifier": "Modifier",
            "rescue_priority_score": "Priority Score",
        }
        table_df = table_df.rename(columns=rename_map)
        if "Modifier" in table_df.columns:
            table_df["Modifier"] = table_df["Modifier"].map(lambda x: f"{x:.2f}")
        if "Priority Score" in table_df.columns:
            table_df["Priority Score"] = table_df["Priority Score"].map(lambda x: f"{x:.3f}")
        table = table_ax.table(
            cellText=table_df.values,
            colLabels=table_df.columns,
            loc="center",
            cellLoc="left",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.3)
    table_ax.set_title("Top-10 Tractability Snapshot")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig6_benchmarking.png", dpi=150)
    plt.savefig(FIGURES_DIR / "fig6_benchmarking.svg")
    plt.close(fig)


def save_alphamissense_concordance_figure(
    concordance_df: pd.DataFrame,
    concordance_summary: pd.DataFrame,
    auc_df: pd.DataFrame,
) -> None:
    if concordance_df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    scatter_ax = axes[0]
    for group, color in CLINVAR_GROUP_COLORS.items():
        sub = concordance_df[concordance_df["clinvar_group"] == group]
        if sub.empty:
            continue
        scatter_ax.scatter(
            sub["alphamissense_score"],
            sub["structural_opportunity_score"],
            s=18,
            alpha=0.55,
            c=color,
            label=group,
        )
    if not concordance_summary.empty:
        row = concordance_summary.iloc[0]
        scatter_ax.axvline(row["alphamissense_high_cutoff"], color="#6C757D", linestyle="--", linewidth=1)
        scatter_ax.axhline(row["structural_high_cutoff"], color="#6C757D", linestyle="--", linewidth=1)
        scatter_ax.text(
            0.03,
            0.97,
            f"Spearman rho={row['spearman_rho']:.2f}\nPearson r={row['pearson_r']:.2f}\nN={int(row['n_variants'])}",
            transform=scatter_ax.transAxes,
            va="top",
        )
    scatter_ax.set_title("AlphaMissense vs Structural Opportunity")
    scatter_ax.set_xlabel("AlphaMissense Score")
    scatter_ax.set_ylabel("Structural Opportunity Score")
    scatter_ax.legend(frameon=False, fontsize=8, loc="lower right")

    bar_ax = axes[1]
    use = auc_df[auc_df["endpoint"].isin(["responsive_any", "responsive_current50", "reaches_wt_current"])].copy()
    if not use.empty:
        order = ["alphamissense_only", "rescue_priority_only", "combined_mean_z"]
        endpoint_order = ["responsive_any", "responsive_current50", "reaches_wt_current"]
        use["model_name"] = pd.Categorical(use["model_name"], categories=order, ordered=True)
        use["endpoint"] = pd.Categorical(use["endpoint"], categories=endpoint_order, ordered=True)
        sns.barplot(data=use, x="endpoint", y="auc", hue="model_name", ax=bar_ax)
        bar_ax.set_ylim(0, 1.0)
        bar_ax.set_xlabel("")
        bar_ax.set_ylabel("Raw Score AUC")
        bar_ax.set_title("Bhatt Retigabine Benchmark")
        bar_ax.set_xticklabels(["Responsive\n(any)", "Responsive\n(>=50% current)", "Strict WT\ncurrent"])
        handles, _ = bar_ax.get_legend_handles_labels()
        bar_ax.legend(handles, ["AlphaMissense", "Rescue priority", "Combined z-mean"], frameon=False, fontsize=8, title="")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "supp_alphamissense_concordance.png", dpi=150)
    plt.savefig(FIGURES_DIR / "supp_alphamissense_concordance.svg")
    plt.close(fig)


def save_operating_point_confusion_figure(metrics_df: pd.DataFrame, confusion_df: pd.DataFrame) -> None:
    if metrics_df.empty or confusion_df.empty:
        return
    models = ["rescue_priority_score", "alphamissense_score"]
    available = [m for m in models if m in confusion_df["model_name"].unique()]
    if not available:
        return
    fig, axes = plt.subplots(1, len(available), figsize=(6 * len(available), 5))
    if len(available) == 1:
        axes = [axes]
    for ax, model_name in zip(axes, available):
        sub = confusion_df[confusion_df["model_name"] == model_name].copy()
        pivot = (
            sub.pivot(index="actual_class", columns="predicted_class", values="count")
            .reindex(index=["positive", "negative"], columns=["positive", "negative"])
            .fillna(0)
        )
        sns.heatmap(pivot, annot=True, fmt=".0f", cmap="Blues", cbar=False, ax=ax, square=True)
        metric_row = metrics_df[metrics_df["model_name"] == model_name]
        title = model_name.replace("_", " ")
        if not metric_row.empty:
            row = metric_row.iloc[0]
            title = f"{title}\nMCC={row['mcc']:.2f}, F1={row['f1_score']:.2f}"
        ax.set_title(title)
        ax.set_xlabel("Predicted class")
        ax.set_ylabel("Actual class")
        ax.set_xticklabels(["Positive", "Negative"], rotation=0)
        ax.set_yticklabels(["Positive", "Negative"], rotation=0)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "supp_operating_point_confusion.png", dpi=150)
    plt.savefig(FIGURES_DIR / "supp_operating_point_confusion.svg")
    plt.close(fig)


def save_modular_benchmark_figure(
    kcnq2_overlap: pd.DataFrame,
    modular_auc_df: pd.DataFrame,
) -> None:
    if kcnq2_overlap.empty:
        return

    use = kcnq2_overlap.copy()
    if "revel_score" not in use.columns:
        revel_path = ROOT / "results_final" / "exports" / "kcnq2_bhatt_revel_scores.csv"
        if revel_path.exists():
            revel = pd.read_csv(revel_path)[["protein_change", "revel_score"]]
            use = use.merge(revel, on="protein_change", how="left")
    use["am_x_tractability"] = use["alphamissense_score"] * use["tractability_modifier"]
    if "revel_score" in use.columns:
        use["revel_x_tractability"] = use["revel_score"] * use["tractability_modifier"]

    plot_models = [
        ("revel_x_tractability", "REVEL × tractability"),
        ("revel_only", "REVEL"),
        ("am_x_tractability", "AlphaMissense × tractability"),
        ("alphamissense_only", "AlphaMissense"),
        ("structural_only", "Structural opportunity"),
        ("rescue_priority_only", "Rescue priority"),
    ]
    label_map = dict(plot_models)
    score_map = {
        "revel_x_tractability": "revel_x_tractability",
        "revel_only": "revel_score",
        "am_x_tractability": "am_x_tractability",
        "alphamissense_only": "alphamissense_score",
        "structural_only": "structural_opportunity_score",
        "rescue_priority_only": "rescue_priority_score",
    }

    auc_rows = []
    for endpoint_name, endpoint_col in [
        ("reaches_wt_current", "bhatt_retig_reaches_wt_current"),
        ("responsive_current50", "bhatt_retig_responsive_current50"),
        ("responsive_any", "bhatt_retig_responsive_any"),
    ]:
        for model_name, score_col in score_map.items():
            if score_col not in use.columns:
                continue
            sub = use.dropna(subset=[score_col, endpoint_col]).copy()
            if len(sub) < 3 or len(sub[endpoint_col].astype(int).unique()) < 2:
                continue
            fpr, tpr, _ = roc_curve(sub[endpoint_col].astype(int), sub[score_col].astype(float))
            auc_rows.append(
                {
                    "endpoint": endpoint_name,
                    "model_name": model_name,
                    "score_column": score_col,
                    "n_variants_used": int(len(sub)),
                    "auc": float(auc(fpr, tpr)),
                }
            )
    auc_df = pd.DataFrame(auc_rows)
    if not modular_auc_df.empty:
        modular_auc_df = modular_auc_df.copy()

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw={"width_ratios": [1, 1, 1.05]})

    roc_specs = [
        ("bhatt_retig_reaches_wt_current", "Strict WT-current Endpoint"),
        ("bhatt_retig_responsive_any", "Broad Responsive Endpoint"),
    ]
    for ax, (endpoint_col, title) in zip(axes[:2], roc_specs):
        for model_name, legend_label in plot_models:
            score_col = score_map.get(model_name)
            if score_col not in use.columns:
                continue
            roc_df = use.dropna(subset=[endpoint_col, score_col]).copy()
            if len(roc_df) < 3 or len(roc_df[endpoint_col].astype(int).unique()) < 2:
                continue
            fpr, tpr, _ = roc_curve(roc_df[endpoint_col].astype(int), roc_df[score_col].astype(float))
            auc_val = auc(fpr, tpr)
            ax.plot(
                fpr,
                tpr,
                color=MODEL_COLORS.get(model_name, "#6C757D"),
                linewidth=2.5,
                label=f"{legend_label} ({auc_val:.3f})",
            )
        ax.plot([0, 1], [0, 1], color="#D3D3D3", linestyle="--", linewidth=1)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(title)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.legend(frameon=False, fontsize=8, loc="lower right")

    bar_ax = axes[2]
    endpoint_order = ["reaches_wt_current", "responsive_current50", "responsive_any"]
    model_order = [m for m, _ in plot_models if m in auc_df["model_name"].unique()]
    plot_df = auc_df[auc_df["model_name"].isin(model_order)].copy()
    plot_df["endpoint"] = pd.Categorical(plot_df["endpoint"], categories=endpoint_order, ordered=True)
    plot_df["model_name"] = pd.Categorical(plot_df["model_name"], categories=model_order, ordered=True)
    sns.barplot(
        data=plot_df,
        x="endpoint",
        y="auc",
        hue="model_name",
        palette={m: MODEL_COLORS.get(m, "#6C757D") for m in model_order},
        ax=bar_ax,
    )
    bar_ax.set_ylim(0, 1.02)
    bar_ax.set_xlabel("")
    bar_ax.set_ylabel("ROC-AUC")
    bar_ax.set_title("AUC Summary Across Rescue Endpoints")
    bar_ax.set_xticklabels(["Strict WT\ncurrent", ">=50% current", "Responsive\n(any)"])
    handles, labels = bar_ax.get_legend_handles_labels()
    bar_ax.legend(
        handles,
        [label_map.get(lbl, lbl.replace("_", " ")) for lbl in labels],
        frameon=False,
        fontsize=8,
        title="",
        loc="lower right",
    )

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "supp_modular_benchmark.png", dpi=150)
    plt.savefig(FIGURES_DIR / "supp_modular_benchmark.svg")
    plt.close(fig)


def _parse_first_model_pdbqt(pdbqt_path: str) -> tuple[list[dict], float | None]:
    atoms = []
    best_dg = None
    in_model = False
    with open(pdbqt_path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if in_model:
                    break
                in_model = line.strip().split()[-1] == "1"
                continue
            if not in_model:
                continue
            if line.startswith("REMARK VINA RESULT:") and best_dg is None:
                try:
                    best_dg = float(line.split()[3])
                except Exception:
                    best_dg = None
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[77:79].strip() or atom_name[0]
                atoms.append({"atom_name": atom_name, "x": x, "y": y, "z": z, "element": element.title()})
    return atoms, best_dg


def _infer_bonds(atoms: list[dict]) -> list[tuple[int, int]]:
    radii = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05, "P": 1.07}
    bonds = []
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=float)
    if len(coords) < 2:
        return bonds
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            ei = atoms[i]["element"]
            ej = atoms[j]["element"]
            cutoff = radii.get(ei, 0.8) + radii.get(ej, 0.8) + 0.45
            d = float(np.linalg.norm(coords[i] - coords[j]))
            if 0.4 < d <= cutoff:
                bonds.append((i, j))
    return bonds


def _atom_color(element: str) -> str:
    return {
        "C": "#2F4858",
        "N": "#3D5A80",
        "O": "#E63946",
        "S": "#E9C46A",
        "F": "#2A9D8F",
        "P": "#9D4EDD",
        "H": "#BDBDBD",
    }.get(element, "#6C757D")


def _plot_atom_sticks(ax, atoms: list[dict], bonds: list[tuple[int, int]], atom_scale: int = 36, override_color: str | None = None) -> None:
    for i, j in bonds:
        xyz = np.array(
            [
                [atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]],
                [atoms[j]["x"], atoms[j]["y"], atoms[j]["z"]],
            ]
        )
        line_color = override_color or "#5C677D"
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color=line_color, linewidth=2.0, alpha=0.95)
    for atom in atoms:
        color = override_color or _atom_color(atom["element"])
        ax.scatter(atom["x"], atom["y"], atom["z"], s=atom_scale, color=color, depthshade=False, edgecolor="none")


def _load_best_chain_atoms(structure_path: str, residue_num: int, ligand_center: np.ndarray, radius: float = 10.0) -> tuple[list[dict], list[dict]]:
    structure = PDBParser(QUIET=True).get_structure("rec", structure_path)
    model = structure[0]
    chains = {c.id: [r for r in c if r.get_id()[0] == " "] for c in model}
    best_chain_id = max(chains, key=lambda cid: len(chains[cid]))
    pocket_atoms = []
    residue_atoms = []
    for residue in chains[best_chain_id]:
        rid = int(residue.id[1])
        for atom in residue:
            if atom.element == "H":
                continue
            coord = atom.coord.astype(float)
            atom_dict = {
                "atom_name": atom.get_name().strip(),
                "x": float(coord[0]),
                "y": float(coord[1]),
                "z": float(coord[2]),
                "element": (atom.element or atom.get_name()[0]).title(),
            }
            if np.linalg.norm(coord - ligand_center) <= radius:
                pocket_atoms.append(atom_dict)
            if rid == residue_num:
                residue_atoms.append(atom_dict)
    return pocket_atoms, residue_atoms


def save_pose_render_figure() -> None:
    specs = [
        {
            "label": "KCNQ2 T276I",
            "drug": "Retigabine",
            "structure_path": ROOT / "pipeline_data/structures/KCNQ2_7CR0.pdb",
            "pose_path": ROOT / "pipeline_data/docking/top10/KCNQ2_T276I_Retigabine_Ezogabine_.pdbqt",
            "residue_num": 276,
        },
        {
            "label": "KCNQ1 L273I",
            "drug": "ML277",
            "structure_path": ROOT / "pipeline_data/structures/KCNQ1_7XNN.pdb",
            "pose_path": ROOT / "pipeline_data/docking/matrix/KCNQ1_L273I_ML277.pdbqt",
            "residue_num": 273,
        },
        {
            "label": "KCNQ4 L281S",
            "drug": "Retigabine",
            "structure_path": ROOT / "pipeline_data/structures/KCNQ4_7BYL.pdb",
            "pose_path": ROOT / "pipeline_data/docking/top10/KCNQ4_L281S_Retigabine_Ezogabine_.pdbqt",
            "residue_num": 281,
        },
    ]
    fig = plt.figure(figsize=(15, 5))
    for idx, spec in enumerate(specs, start=1):
        if not spec["structure_path"].exists() or not spec["pose_path"].exists():
            continue
        ligand_atoms, best_dg = _parse_first_model_pdbqt(str(spec["pose_path"]))
        if not ligand_atoms:
            continue
        ligand_center = np.mean(np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms]), axis=0)
        pocket_atoms, residue_atoms = _load_best_chain_atoms(str(spec["structure_path"]), spec["residue_num"], ligand_center)
        ligand_bonds = _infer_bonds(ligand_atoms)
        residue_bonds = _infer_bonds(residue_atoms)

        ax = fig.add_subplot(1, 3, idx, projection="3d")
        if pocket_atoms:
            pocket_xyz = np.array([[a["x"], a["y"], a["z"]] for a in pocket_atoms])
            ax.scatter(
                pocket_xyz[:, 0],
                pocket_xyz[:, 1],
                pocket_xyz[:, 2],
                s=10,
                color="#C9D1D9",
                alpha=0.15,
                depthshade=False,
                edgecolor="none",
            )
        _plot_atom_sticks(ax, ligand_atoms, ligand_bonds, atom_scale=48)
        _plot_atom_sticks(ax, residue_atoms, residue_bonds, atom_scale=58, override_color="#F77F00")

        all_xyz = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms + residue_atoms + pocket_atoms[:150]])
        center = all_xyz.mean(axis=0)
        span = max(8.0, float(np.max(np.ptp(all_xyz, axis=0))) / 2 + 2.0)
        ax.set_xlim(center[0] - span, center[0] + span)
        ax.set_ylim(center[1] - span, center[1] + span)
        ax.set_zlim(center[2] - span, center[2] + span)
        ax.view_init(elev=18, azim=-55)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_box_aspect((1, 1, 1))
        ax.set_title(f"{spec['label']}\n{spec['drug']} | ΔG {best_dg:.2f} kcal/mol" if best_dg is not None else f"{spec['label']}\n{spec['drug']}")
        ax.text2D(0.03, 0.04, f"Residue {spec['residue_num']}", transform=ax.transAxes, color="#B34700", fontsize=9)

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#C9D1D9", alpha=0.45, markersize=8, label="Pocket atoms"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2F4858", markersize=8, label="Ligand"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#F77F00", markersize=8, label="Candidate residue"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3, frameon=False)
    plt.tight_layout(rect=(0, 0.05, 1, 1))
    plt.savefig(FIGURES_DIR / "fig7_pose_renders.png", dpi=220)
    plt.savefig(FIGURES_DIR / "fig7_pose_renders.svg")
    plt.close(fig)


def save_sensitivity_heatmap(sensitivity_effects: pd.DataFrame) -> None:
    if sensitivity_effects.empty:
        return
    pivot = sensitivity_effects.pivot(index="candidate", columns="weight_name", values="median_rank_span")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(pivot))))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="YlOrRd", ax=ax)
    ax.set_title("Candidate Rank Sensitivity by Weight")
    ax.set_xlabel("Weight Component")
    ax.set_ylabel("Candidate")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "supp_weight_sensitivity.png", dpi=150)
    plt.savefig(FIGURES_DIR / "supp_weight_sensitivity.svg")
    plt.close(fig)
