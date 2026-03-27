from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from kcnq_pipeline.config import EXPORT_DIR


def _extract_nested(obj: dict, path: str):
    cur = obj
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def _pick_revel_score(hit: dict) -> tuple[float | None, float | None]:
    candidate_paths = [
        "dbnsfp.revel.score",
        "dbnsfp.revel.rankscore",
        "dbnsfp.REVEL.score",
        "dbnsfp.REVEL.rankscore",
    ]
    found = {}
    for path in candidate_paths:
        value = _extract_nested(hit, path)
        if isinstance(value, list):
            value = value[0] if value else None
        if value is not None:
            found[path] = value
    score = found.get("dbnsfp.revel.score") or found.get("dbnsfp.REVEL.score")
    rankscore = found.get("dbnsfp.revel.rankscore") or found.get("dbnsfp.REVEL.rankscore")
    try:
        score = float(score) if score is not None else None
    except Exception:
        score = None
    try:
        rankscore = float(rankscore) if rankscore is not None else None
    except Exception:
        rankscore = None
    return score, rankscore


def main() -> None:
    overlap = pd.read_csv(EXPORT_DIR / "kcnq2_bhatt_retigabine_overlap.csv")
    overlap = overlap.drop_duplicates("protein_change").copy()
    overlap["hgvs_g"] = overlap.apply(
        lambda r: f"chr{str(r['chromosome']).replace('chr', '')}:g.{int(r['position'])}{r['ref']}>{r['alt']}",
        axis=1,
    )

    session = requests.Session()
    session.headers.update({"User-Agent": "codex-kcnq-revel-fetch/1.0"})
    rows = []
    for _, row in overlap.iterrows():
        hgvs_g = row["hgvs_g"]
        url = f"https://myvariant.info/v1/variant/{hgvs_g}"
        resp = session.get(url, params={"fields": "dbnsfp.revel", "assembly": "hg38"}, timeout=60)
        payload = resp.json()
        score, rankscore = _pick_revel_score(payload if isinstance(payload, dict) else {})
        rows.append(
            {
                "protein_change": row["protein_change"],
                "hgvs_g": hgvs_g,
                "revel_score": score,
                "revel_rankscore": rankscore,
                "found": score is not None or rankscore is not None,
                "raw_response_snippet": json.dumps(payload)[:500],
            }
        )
        time.sleep(0.1)

    out = pd.DataFrame(rows)
    out.to_csv(EXPORT_DIR / "kcnq2_bhatt_revel_scores.csv", index=False)
    print(out[["found", "revel_score", "revel_rankscore"]].notna().sum().to_dict())
    print(EXPORT_DIR / "kcnq2_bhatt_revel_scores.csv")


if __name__ == "__main__":
    main()
