from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select

from .config import DOMAIN_BOUNDARIES, REGION_LABELS, SELECTED_STRUCTURES, STRUCTURE_DIR, UNIPROT_IDS
from .utils import download, wait_between_requests


def fetch_structures() -> dict[str, dict]:
    files: dict[str, dict] = {}
    for gene, meta in SELECTED_STRUCTURES.items():
        if meta["source"] == "PDB":
            pdb_id = meta["pdb_id"]
            dest = STRUCTURE_DIR / f"{gene}_{pdb_id}.pdb"
            if not dest.exists():
                try:
                    download(f"https://files.rcsb.org/download/{pdb_id}.pdb", dest)
                except Exception:
                    cif_dest = STRUCTURE_DIR / f"{gene}_{pdb_id}.cif"
                    download(f"https://files.rcsb.org/download/{pdb_id}.cif", cif_dest)
                    dest = cif_dest
            files[gene] = {**meta, "path": str(dest)}
        else:
            uniprot = meta["uniprot"]
            dest = STRUCTURE_DIR / f"{gene}_AF_{uniprot}.pdb"
            if not dest.exists():
                api = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}", timeout=60)
                api.raise_for_status()
                pdb_url = api.json()[0]["pdbUrl"]
                download(pdb_url, dest)
            files[gene] = {**meta, "path": str(dest)}
        wait_between_requests(0.15)
    return files


def _load_structure(path: str, name: str):
    p = Path(path)
    if p.suffix.lower() == ".cif":
        return MMCIFParser(QUIET=True).get_structure(name, str(p))
    return PDBParser(QUIET=True).get_structure(name, str(p))


def parse_structure_coords(structure_files: dict[str, dict]) -> tuple[dict[str, dict], dict[str, tuple[float, float, float] | None]]:
    structure_coords: dict[str, dict] = {}
    pocket_centroids: dict[str, tuple[float, float, float] | None] = {}
    for gene, meta in structure_files.items():
        structure = _load_structure(meta["path"], gene)
        model = structure[0]
        chains = {c.id: [r for r in c if r.get_id()[0] == " "] for c in model}
        best_chain = max(chains, key=lambda cid: len(chains[cid]))
        coords = {}
        for residue in chains[best_chain]:
            if "CA" in residue:
                coords[int(residue.id[1])] = tuple(float(x) for x in residue["CA"].coord)
        structure_coords[gene] = {"chain": best_chain, "coords": coords}
        lo, hi = DOMAIN_BOUNDARIES[gene]["Ligand_pocket"]
        points = [np.array(coords[r]) for r in range(lo, hi + 1) if r in coords]
        pocket_centroids[gene] = tuple(np.mean(points, axis=0)) if points else None
    return structure_coords, pocket_centroids


def assign_region(gene: str, residue_num: float | int | None) -> str:
    if pd.isna(residue_num):
        return "Unknown"
    residue_num = int(residue_num)
    for region, (lo, hi) in DOMAIN_BOUNDARIES[gene].items():
        if lo <= residue_num <= hi:
            return region
    return "Unknown"


def dist_to_centroid(coord, centroid) -> float | None:
    if coord is None or centroid is None:
        return np.nan
    return float(np.linalg.norm(np.array(coord) - np.array(centroid)))


def annotate_variants(df: pd.DataFrame, structure_coords: dict[str, dict], pocket_centroids: dict[str, tuple[float, float, float] | None]) -> pd.DataFrame:
    out = df.copy()
    out["structural_region"] = out.apply(lambda r: assign_region(r["gene"], r["residue_num"]), axis=1)
    out["region_label"] = out["structural_region"].map(REGION_LABELS)

    def get_residue_coord(gene: str, residue_num):
        if pd.isna(residue_num):
            return None
        coords = structure_coords.get(gene, {}).get("coords", {})
        residue_num = int(residue_num)
        for offset in [0, 1, -1, 2, -2]:
            if residue_num + offset in coords:
                return coords[residue_num + offset]
        return None

    out["pocket_dist_A"] = out.apply(
        lambda r: dist_to_centroid(get_residue_coord(r["gene"], r["residue_num"]), pocket_centroids.get(r["gene"])),
        axis=1,
    )
    out["structural_coords_available"] = out["pocket_dist_A"].notna()
    return out


class PocketSelect(Select):
    def __init__(self, res_range: tuple[int, int], chain_id: str):
        self.lo, self.hi = res_range
        self.chain_id = chain_id

    def accept_residue(self, residue):
        rid = residue.get_id()[1]
        return (self.lo - 30) <= rid <= (self.hi + 30)

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def export_pocket_pdb(gene: str, structure_files: dict[str, dict], structure_coords: dict[str, dict], out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = structure_files[gene]
    structure = _load_structure(meta["path"], gene)
    chain = structure_coords[gene]["chain"]
    out_pdb = out_dir / f"{gene}_pocket.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb), PocketSelect(DOMAIN_BOUNDARIES[gene]["Ligand_pocket"], chain))
    return out_pdb

