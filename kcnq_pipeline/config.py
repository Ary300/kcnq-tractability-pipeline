from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "pipeline_data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
BENCHMARK_DIR = DATA_DIR / "benchmarks"
STRUCTURE_DIR = DATA_DIR / "structures"
DOCKING_DIR = DATA_DIR / "docking"
RESULTS_DIR = ROOT / "results_final"
FIGURES_DIR = RESULTS_DIR / "figures"
EXPORT_DIR = RESULTS_DIR / "exports"
CACHE_DIR = DATA_DIR / "cache"

for directory in [
    DATA_DIR,
    RAW_DIR,
    INTERIM_DIR,
    BENCHMARK_DIR,
    STRUCTURE_DIR,
    DOCKING_DIR,
    RESULTS_DIR,
    FIGURES_DIR,
    EXPORT_DIR,
    CACHE_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

GENES = ["KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5"]

GENE_IDS = {
    "KCNQ1": "ENSG00000053918",
    "KCNQ2": "ENSG00000075043",
    "KCNQ3": "ENSG00000184156",
    "KCNQ4": "ENSG00000117601",
    "KCNQ5": "ENSG00000185760",
}

UNIPROT_IDS = {
    "KCNQ1": "P51787",
    "KCNQ2": "O43526",
    "KCNQ3": "O43525",
    "KCNQ4": "P56696",
    "KCNQ5": "Q9NR82",
}

SELECTED_STRUCTURES = {
    "KCNQ1": {"source": "PDB", "pdb_id": "7XNN", "resolution_A": 2.5},
    "KCNQ2": {"source": "PDB", "pdb_id": "7CR0", "resolution_A": 3.1},
    "KCNQ3": {"source": "AlphaFold", "uniprot": "O43525", "pdb_id": "AF-O43525-F1"},
    "KCNQ4": {"source": "PDB", "pdb_id": "7BYL", "resolution_A": 2.5},
    "KCNQ5": {"source": "PDB", "pdb_id": "9J38", "resolution_A": 2.4},
}

DOMAIN_BOUNDARIES = {
    "KCNQ1": {
        "Unknown": (1, 59),
        "VSD": (60, 200),
        "Pore": (201, 330),
        "Selectivity": (270, 285),
        "Ligand_pocket": (250, 340),
        "CaM": (331, 500),
        "Assembly": (501, 700),
    },
    "KCNQ2": {
        "Unknown": (1, 59),
        "VSD": (60, 200),
        "Pore": (201, 330),
        "Selectivity": (270, 285),
        "Ligand_pocket": (201, 330),
        "CaM": (331, 500),
        "Assembly": (501, 900),
    },
    "KCNQ3": {
        "Unknown": (1, 59),
        "VSD": (60, 200),
        "Pore": (201, 330),
        "Selectivity": (270, 285),
        "Ligand_pocket": (201, 330),
        "CaM": (331, 500),
        "Assembly": (501, 900),
    },
    "KCNQ4": {
        "Unknown": (1, 59),
        "VSD": (60, 200),
        "Pore": (201, 330),
        "Selectivity": (270, 285),
        "Ligand_pocket": (201, 330),
        "CaM": (331, 500),
        "Assembly": (501, 700),
    },
    "KCNQ5": {
        "Unknown": (1, 59),
        "VSD": (60, 200),
        "Pore": (201, 330),
        "Selectivity": (270, 285),
        "Ligand_pocket": (201, 330),
        "CaM": (331, 500),
        "Assembly": (501, 700),
    },
}

REGION_LABELS = {
    "Unknown": "Unknown/Unresolved",
    "VSD": "Voltage-Sensing Domain",
    "Pore": "Pore Domain (S5-S6)",
    "Selectivity": "Pore/Selectivity Filter",
    "Ligand_pocket": "Ligand-Binding Pocket",
    "CaM": "CaM/Interface Region",
    "Assembly": "Trafficking/Assembly",
}

KNOWN_KCNQ_DRUGS = [
    {
        "gene": "KCNQ1",
        "drug_name": "Mexiletine",
        "max_phase": 4,
        "binding_region": "Pore Domain (S5-S6)",
        "smiles": "CC(Nc1c(C)cccc1C)CN",
    },
    {
        "gene": "KCNQ1",
        "drug_name": "Flecainide",
        "max_phase": 4,
        "binding_region": "Pore Domain (S5-S6)",
        "smiles": "OCC(F)(F)OC1=CC(=CC(=C1)C(=O)NCC2CCCCN2)OCC(F)(F)O",
    },
    {
        "gene": "KCNQ1",
        "drug_name": "ML277",
        "max_phase": 1,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "O=C(Nc1ccc(S(=O)(=O)N2CCOCC2)cc1)c1ccc(F)cc1",
    },
    {
        "gene": "KCNQ2",
        "drug_name": "Retigabine (Ezogabine)",
        "max_phase": 4,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CCOC(=O)Nc1ccc(NCC2=CC=CC=C2F)cc1N",
    },
    {
        "gene": "KCNQ2",
        "drug_name": "HN37",
        "max_phase": 2,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CC1=CC(=CC(=C1NC(=O)OC)C)N(CC#C)CC2=CC=C(C=C2)F",
    },
    {
        "gene": "KCNQ2",
        "drug_name": "XEN1101",
        "max_phase": 3,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CC1=CC(=CC(=C1NC(=O)CC(C)(C)C)C)N2CCC3=C(C2)C=CC(=C3)F",
    },
    {
        "gene": "KCNQ3",
        "drug_name": "Retigabine (Ezogabine)",
        "max_phase": 4,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CCOC(=O)Nc1ccc(NCC2=CC=CC=C2F)cc1N",
    },
    {
        "gene": "KCNQ3",
        "drug_name": "XEN1101",
        "max_phase": 3,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CC1=CC(=CC(=C1NC(=O)CC(C)(C)C)C)N2CCC3=C(C2)C=CC(=C3)F",
    },
    {
        "gene": "KCNQ4",
        "drug_name": "Retigabine (Ezogabine)",
        "max_phase": 4,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CCOC(=O)Nc1ccc(NCC2=CC=CC=C2F)cc1N",
    },
    {
        "gene": "KCNQ5",
        "drug_name": "Retigabine (Ezogabine)",
        "max_phase": 4,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CCOC(=O)Nc1ccc(NCC2=CC=CC=C2F)cc1N",
    },
    {
        "gene": "KCNQ5",
        "drug_name": "HN37",
        "max_phase": 2,
        "binding_region": "Ligand-Binding Pocket",
        "smiles": "CC1=CC(=CC(=C1NC(=O)OC)C)N(CC#C)CC2=CC=C(C=C2)F",
    },
]

REGION_DRUG_MAP = {
    "Ligand-Binding Pocket": 1.0,
    "Pore/Selectivity Filter": 0.9,
    "Pore Domain (S5-S6)": 0.9,
    "Voltage-Sensing Domain": 0.6,
    "CaM/Interface Region": 0.5,
    "Trafficking/Assembly": 0.3,
    "Unknown/Unresolved": 0.2,
}

CLINVAR_TIER_MAP = {
    "Pathogenic": 1.0,
    "Pathogenic/Likely pathogenic": 0.95,
    "Likely pathogenic": 0.8,
    "Conflicting classifications of pathogenicity": 0.5,
    "Uncertain significance": 0.3,
    "gnomAD_population": 0.1,
    "Likely benign": 0.05,
    "Benign/Likely benign": 0.02,
    "Benign": 0.01,
    "not provided": 0.2,
    "": 0.2,
}

EXCLUDE_SIG = {
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
    "gnomAD_population",
    "not provided",
    "",
}

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
ALPHAMISSENSE_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
CADD_URL = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7"

VINA_CANDIDATES = [
    ROOT.parent / "vina_1.2.6_mac_aarch64",
    ROOT.parent / "vina_1.2.6_mac_x86_64",
]
