from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from .config import DOCKING_DIR, KNOWN_KCNQ_DRUGS, VINA_CANDIDATES
from .structures import export_pocket_pdb


def resolve_vina() -> str:
    found = shutil.which("vina")
    if found:
        return found
    for candidate in VINA_CANDIDATES:
        if candidate.exists():
            candidate.chmod(candidate.stat().st_mode | 0o111)
            return str(candidate)
    raise FileNotFoundError("No Vina binary available")


def write_smiles_table() -> pd.DataFrame:
    df = pd.DataFrame(KNOWN_KCNQ_DRUGS).drop_duplicates(["drug_name"])
    return df[df["smiles"].astype(str).str.len() > 0].copy()


def prepare_ligands(ligands_dir: Path) -> pd.DataFrame:
    ligands_dir.mkdir(parents=True, exist_ok=True)
    df = write_smiles_table()
    prepared = []
    for row in df.itertuples(index=False):
        safe = re.sub(r"[^A-Za-z0-9]+", "_", row.drug_name).strip("_")
        sdf_path = ligands_dir / f"{safe}.sdf"
        pdbqt_path = ligands_dir / f"{safe}.pdbqt"
        mol = Chem.MolFromSmiles(row.smiles)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            continue
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        mol.SetProp("_Name", safe)
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()
        cmd = [
            "mk_prepare_ligand.py",
            "-i",
            str(sdf_path),
            "-o",
            str(pdbqt_path),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        prepared.append({"drug_name": row.drug_name, "pdbqt_path": str(pdbqt_path), "smiles": row.smiles})
    return pd.DataFrame(prepared)


def prepare_receptors(structure_files: dict[str, dict], structure_coords: dict[str, dict], pocket_centroids: dict[str, tuple[float, float, float] | None], receptors_dir: Path) -> tuple[dict[str, str], dict[str, dict]]:
    receptors_dir.mkdir(parents=True, exist_ok=True)
    receptor_paths: dict[str, str] = {}
    boxes: dict[str, dict] = {}
    for gene, centroid in pocket_centroids.items():
        if centroid is None:
            continue
        pocket_pdb = export_pocket_pdb(gene, structure_files, structure_coords, receptors_dir)
        receptor_pdbqt = receptors_dir / f"{gene}_pocket.pdbqt"
        cmd = [
            "obabel",
            "-ipdb",
            str(pocket_pdb),
            "-opdbqt",
            "-xr",
            "-O",
            str(receptor_pdbqt),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        receptor_paths[gene] = str(receptor_pdbqt)
        boxes[gene] = {
            "center_x": float(centroid[0]),
            "center_y": float(centroid[1]),
            "center_z": float(centroid[2]),
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
        }
    return receptor_paths, boxes


def run_vina(vina_bin: str, receptor: str, ligand: str, box: dict, out_path: Path, exhaustiveness: int = 16) -> float | None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        vina_bin,
        "--receptor",
        receptor,
        "--ligand",
        ligand,
        "--center_x",
        str(box["center_x"]),
        "--center_y",
        str(box["center_y"]),
        "--center_z",
        str(box["center_z"]),
        "--size_x",
        str(box["size_x"]),
        "--size_y",
        str(box["size_y"]),
        "--size_z",
        str(box["size_z"]),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        "5",
        "--cpu",
        "4",
        "--out",
        str(out_path),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    for line in proc.stdout.splitlines():
        match = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if match:
            return float(match.group(1))
    return None


def docking_score(dg: float | None) -> float:
    if dg is None or pd.isna(dg):
        return np.nan
    return max(0.0, min(1.0, (float(dg) - (-4.0)) / (-8.0 - (-4.0))))


def run_top10_docking(df_top10: pd.DataFrame, receptor_paths: dict[str, str], boxes: dict[str, dict], ligands_df: pd.DataFrame) -> pd.DataFrame:
    vina_bin = resolve_vina()
    ligand_map = dict(zip(ligands_df["drug_name"], ligands_df["pdbqt_path"]))
    rows = []
    for row in df_top10.itertuples(index=False):
        receptor = receptor_paths.get(row.gene)
        ligand = ligand_map.get(row.candidate_drug)
        box = boxes.get(row.gene)
        dg = np.nan
        note = "ok"
        if receptor and ligand and box:
            out_path = DOCKING_DIR / "top10" / f"{row.gene}_{row.protein_change}_{re.sub(r'[^A-Za-z0-9]+','_',row.candidate_drug)}.pdbqt"
            result = run_vina(vina_bin, receptor, ligand, box, out_path, exhaustiveness=16)
            dg = result if result is not None else np.nan
            if result is None:
                note = "failed"
        else:
            note = "missing_inputs"
        rows.append(
            {
                "gene": row.gene,
                "protein_change": row.protein_change,
                "candidate_drug": row.candidate_drug,
                "dG_WT": dg,
                "docking_wt_score": docking_score(dg),
                "note": note,
            }
        )
    return pd.DataFrame(rows)


def applicable_drugs_for_gene(gene: str) -> list[str]:
    gene_drugs = []
    for row in KNOWN_KCNQ_DRUGS:
        if row["gene"] == gene and row["smiles"]:
            gene_drugs.append(row["drug_name"])
    return sorted(set(gene_drugs))


def run_full_matrix_docking(df_top10: pd.DataFrame, receptor_paths: dict[str, str], boxes: dict[str, dict], ligands_df: pd.DataFrame) -> pd.DataFrame:
    vina_bin = resolve_vina()
    ligand_map = dict(zip(ligands_df["drug_name"], ligands_df["pdbqt_path"]))
    rows = []
    for row in df_top10.itertuples(index=False):
        receptor = receptor_paths.get(row.gene)
        box = boxes.get(row.gene)
        for drug in applicable_drugs_for_gene(row.gene):
            ligand = ligand_map.get(drug)
            dg = np.nan
            note = "ok"
            if receptor and ligand and box:
                out_path = DOCKING_DIR / "matrix" / f"{row.gene}_{row.protein_change}_{re.sub(r'[^A-Za-z0-9]+','_',drug)}.pdbqt"
                result = run_vina(vina_bin, receptor, ligand, box, out_path, exhaustiveness=16)
                dg = result if result is not None else np.nan
                if result is None:
                    note = "failed"
            else:
                note = "missing_inputs"
            rows.append(
                {
                    "gene": row.gene,
                    "protein_change": row.protein_change,
                    "drug_name": drug,
                    "dG_WT": dg,
                    "docking_wt_score": docking_score(dg),
                    "note": note,
                }
            )
    return pd.DataFrame(rows)
