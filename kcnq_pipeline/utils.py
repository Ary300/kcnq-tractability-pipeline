from __future__ import annotations

import csv
import gzip
import io
import json
import math
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "kcnq-rescue-pipeline/1.0"})
    return s


def chunks(items: list[str], size: int) -> Iterable[list[str]]:
    for idx in range(0, len(items), size):
        yield items[idx : idx + size]


def parse_spdi(spdi: str) -> tuple[str | None, int | None, str | None, str | None]:
    if not spdi:
        return None, None, None, None
    try:
        seq_id, zero_based, ref, alt = spdi.split(":")
        chrom = seq_id.split(".")[0].replace("NC_0000", "").replace("NC_000", "")
        chrom = chrom.lstrip("0")
        return chrom, int(zero_based) + 1, ref, alt
    except Exception:
        return None, None, None, None


def save_json(obj: object, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2))


def load_json(path: Path, default: object = None) -> object:
    if not path.exists():
        return default
    return json.loads(path.read_text())


def download(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with session().get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with dest.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    fh.write(chunk)
    return dest


def rarity_score(af: float | None) -> float:
    if af is None or pd.isna(af) or af == 0:
        return 1.0
    return max(0.0, 1.0 + math.log10(af) / 5.0)


def ensure_text_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def read_gz_tsv_subset(path: Path, keys: set[tuple[str, int, str, str]], chunksize: int = 200_000) -> pd.DataFrame:
    keep = []
    usecols = [
        "#CHROM",
        "POS",
        "REF",
        "ALT",
        "transcript_id",
        "protein_variant",
        "am_pathogenicity",
        "am_class",
    ]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=3,
        usecols=usecols,
        chunksize=chunksize,
    ):
        chunk["POS"] = chunk["POS"].astype(int)
        chunk["key"] = list(zip(chunk["#CHROM"].str.replace("chr", "", regex=False), chunk["POS"], chunk["REF"], chunk["ALT"]))
        hit = chunk[chunk["key"].isin(keys)].copy()
        if not hit.empty:
            keep.append(hit)
    if keep:
        out = pd.concat(keep, ignore_index=True)
        return out.drop(columns="key")
    return pd.DataFrame(columns=usecols)


def wait_between_requests(seconds: float = 0.25) -> None:
    time.sleep(seconds)
