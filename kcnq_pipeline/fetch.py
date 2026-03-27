from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .config import CACHE_DIR, ENTREZ_BASE, GENES, GENE_IDS, GNOMAD_URL
from .utils import chunks, parse_spdi, session, wait_between_requests


AA3TO1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}


def parse_protein_change(raw_protein_change: str, variation_name: str) -> tuple[str, str, str, int | None]:
    variation_name = variation_name or ""
    raw_protein_change = raw_protein_change or ""

    hgvs_match = re.search(r"\(p\.([A-Za-z]{3}|\*)(\d+)([A-Za-z]{3}|\*)\)", variation_name)
    if hgvs_match:
        aa_ref = AA3TO1.get(hgvs_match.group(1), hgvs_match.group(1))
        aa_alt = AA3TO1.get(hgvs_match.group(3), hgvs_match.group(3))
        residue_num = int(hgvs_match.group(2))
        return f"{aa_ref}{residue_num}{aa_alt}", aa_ref, aa_alt, residue_num

    protein_changes = [p.strip() for p in raw_protein_change.split(",") if p.strip()]
    for protein_change in protein_changes:
        match = re.match(r"^([A-Za-z*]+)(\d+)([A-Za-z*]+)$", protein_change)
        if match:
            aa_ref = match.group(1)
            residue_num = int(match.group(2))
            aa_alt = match.group(3)
            return protein_change, aa_ref, aa_alt, residue_num
    return "", "", "", None


def fetch_clinvar_variants(retmax: int = 10000) -> pd.DataFrame:
    s = session()
    all_records: list[dict] = []
    for gene in GENES:
        search = s.get(
            f"{ENTREZ_BASE}/esearch.fcgi",
            params={
                "db": "clinvar",
                "term": f'{gene}[gene] AND "missense variant"[molecular consequence]',
                "retmax": retmax,
                "retmode": "json",
            },
            timeout=60,
        )
        search.raise_for_status()
        data = search.json()["esearchresult"]
        ids = data["idlist"]
        for batch in chunks(ids, 100):
            summary = s.get(
                f"{ENTREZ_BASE}/esummary.fcgi",
                params={"db": "clinvar", "id": ",".join(batch), "retmode": "json"},
                timeout=60,
            )
            summary.raise_for_status()
            result = summary.json()["result"]
            for uid in batch:
                doc = result.get(uid, {})
                vset = doc.get("variation_set", [])
                if not vset:
                    continue
                raw_pc = doc.get("protein_change", "")
                variation_name = vset[0].get("variation_name", "")
                protein_change, aa_ref, aa_alt, residue_num = parse_protein_change(raw_pc, variation_name)
                clin_sig = doc.get("germline_classification", {}).get("description", "") or doc.get(
                    "clinical_significance", {}
                ).get("description", "")
                trait_set = doc.get("germline_classification", {}).get("trait_set", [])
                canonical_spdi = vset[0].get("canonical_spdi", "")
                chrom, pos, ref, alt = parse_spdi(canonical_spdi)
                all_records.append(
                    {
                        "gene": gene,
                        "clinvar_id": uid,
                        "variant_id": str(uid),
                        "protein_change": protein_change,
                        "aa_ref": aa_ref,
                        "aa_alt": aa_alt,
                        "residue_num": residue_num,
                        "clinical_significance": clin_sig,
                        "trait": trait_set[0].get("trait_name", "") if trait_set else "",
                        "hgvs_c": variation_name,
                        "canonical_spdi": canonical_spdi,
                        "chromosome": chrom,
                        "position": pos,
                        "ref": ref,
                        "alt": alt,
                        "source": "ClinVar",
                    }
                )
            wait_between_requests(0.35)
    df = pd.DataFrame(all_records)
    df = df[df["protein_change"].str.len() > 0].copy()
    df = df[df["residue_num"].notna()].copy()
    return df.reset_index(drop=True)


def fetch_gnomad_variants() -> pd.DataFrame:
    query = """
    query GnomadVariants($geneId: String!, $dataset: DatasetId!) {
      gene(gene_id: $geneId, reference_genome: GRCh38) {
        variants(dataset: $dataset) {
          variant_id
          pos
          ref
          alt
          consequence
          hgvsc
          hgvsp
          exome { ac an af }
          genome { ac an af }
        }
      }
    }
    """
    s = session()
    records: list[dict] = []
    for gene, gene_id in GENE_IDS.items():
        resp = s.post(
            GNOMAD_URL,
            json={"query": query, "variables": {"geneId": gene_id, "dataset": "gnomad_r4"}},
            timeout=120,
            headers={"Content-Type": "application/json"},
        )
        resp.raise_for_status()
        variants = resp.json().get("data", {}).get("gene", {}).get("variants", [])
        for var in variants:
            if var.get("consequence") != "missense_variant":
                continue
            hgvsp = var.get("hgvsp") or ""
            match = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})", hgvsp)
            aa_ref, aa_alt, residue_num, protein_change = "", "", None, ""
            if match:
                aa_ref = AA3TO1.get(match.group(1), match.group(1))
                aa_alt = AA3TO1.get(match.group(3), match.group(3))
                residue_num = int(match.group(2))
                protein_change = f"{aa_ref}{residue_num}{aa_alt}"
            af = None
            ac = None
            if var.get("exome") and var["exome"].get("af") is not None:
                af = var["exome"]["af"]
                ac = var["exome"]["ac"]
            elif var.get("genome") and var["genome"].get("af") is not None:
                af = var["genome"]["af"]
                ac = var["genome"]["ac"]
            records.append(
                {
                    "gene": gene,
                    "variant_id": var.get("variant_id", ""),
                    "protein_change": protein_change,
                    "aa_ref": aa_ref,
                    "aa_alt": aa_alt,
                    "residue_num": residue_num,
                    "clinical_significance": "gnomAD_population",
                    "trait": "",
                    "hgvs_c": var.get("hgvsc", ""),
                    "chromosome": str(var.get("variant_id", "")).split("-")[0].replace("chr", ""),
                    "position": var.get("pos"),
                    "ref": var.get("ref"),
                    "alt": var.get("alt"),
                    "gnomad_af": af,
                    "gnomad_ac": ac,
                    "source": "gnomAD",
                }
            )
        wait_between_requests(0.75)
    df = pd.DataFrame(records)
    df = df[df["protein_change"].str.len() > 0].copy()
    return df.reset_index(drop=True)


def build_variant_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clinvar = fetch_clinvar_variants()
    gnomad = fetch_gnomad_variants()
    clinvar["gnomad_af"] = pd.NA
    clinvar["gnomad_ac"] = pd.NA
    merged = pd.concat([clinvar, gnomad], ignore_index=True, sort=False)
    merged["_key"] = merged["gene"] + "_" + merged["protein_change"]
    merged = merged.sort_values(["source", "variant_id"]).drop_duplicates("_key", keep="first").drop(columns="_key")
    merged = merged.reset_index(drop=True)
    clinvar_only = merged[(merged["source"] == "ClinVar") & (~merged["clinical_significance"].isin({"", "not provided"}))].copy()
    return clinvar.reset_index(drop=True), gnomad.reset_index(drop=True), merged.reset_index(drop=True)
