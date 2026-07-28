"""Offline, zero-heavy-dependency lookup of pre-scored agents.

The common case in practice is "what does *this named drug* score?", and the answer
is already computed and shipped with the package. Two bundled tables are searched,
richest first:

  1. ``data/panel_scored.csv``          - 209 curated toxicology agents; carries drug
                                          class and a short clinical note.
  2. ``data/approved_drugs_scored.csv`` - 2845 approved drugs; broad coverage plus the
                                          within-set amenability rank used in the paper.

Because these tables already contain every number, ``ile-predict "lurasidone"`` answers
instantly **without importing ADMET-AI / PyTorch**. Only a structure that is not in
either table (a novel molecule, a raw SMILES) needs the full prediction stack.
"""
from __future__ import annotations

import functools
import re
from pathlib import Path

import pandas as pd

_DATA = Path(__file__).resolve().parent.parent / "data"


def _norm(name: str) -> str:
    """Normalise a name for matching: lowercase, drop parentheticals and punctuation."""
    s = str(name).lower().strip()
    s = re.sub(r"\(.*?\)", " ", s)            # "ergocalciferol (d2)" -> "ergocalciferol"
    s = re.sub(r"[^a-z0-9]+", " ", s).strip()
    return s


def _category(prob: float) -> str:
    return "High" if prob >= 70 else ("Medium" if prob >= 40 else "Low")


def _ad_from_mw(mw) -> str:
    if mw is None or pd.isna(mw):
        return "unknown"
    return "domain-edge (prediction unreliable)" if (mw > 600 or mw < 100) else "ok"


def _num(row, col):
    if col not in row or pd.isna(row[col]):
        return None
    return round(float(row[col]), 3)


def _clean(x):
    if x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() == "":
        return None
    return str(x).strip()


@functools.lru_cache(maxsize=1)
def _panel() -> pd.DataFrame:
    d = pd.read_csv(_DATA / "panel_scored.csv")
    d = d.rename(columns={
        "Vd_Lkg": "Vd", "ILE_score": "theory_score", "ILE_prob": "ile_prob",
        "AD": "ad_flag", "klinik_not": "note", "class": "cls",
    })
    d["_key"] = d["name"].map(_norm)
    return d


@functools.lru_cache(maxsize=1)
def _approved() -> pd.DataFrame:
    d = pd.read_csv(_DATA / "approved_drugs_scored.csv")
    d["_key"] = d["name"].map(_norm)
    return d


@functools.lru_cache(maxsize=1)
def _aliases() -> dict:
    """Normalised alias -> canonical bundled name (common abbreviations, brand and
    street names, e.g. thc->dronabinol, asa->aspirin, seroquel->quetiapine)."""
    p = _DATA / "aliases.csv"
    if not p.exists():
        return {}
    d = pd.read_csv(p)
    return {_norm(a): str(c) for a, c in zip(d["alias"], d["canonical"])}


def lookup(name: str) -> dict | None:
    """Return a normalised scored record for ``name``, or ``None`` if not bundled.

    The record has stable keys regardless of which table matched::

        query, matched, source, class, logD, Vd, MW, PPB, theory_score,
        ile_prob, category, ad_flag, rank, n_ranked, note
    """
    key = _norm(name)
    if not key:
        return None

    # resolve a common alias (thc, cbd, asa, brand names...) to its canonical name,
    # unless the typed token is itself a bundled entry
    alias_of = None
    if key not in set(_panel()["_key"]) and key not in set(_approved()["_key"]):
        canon = _aliases().get(key)
        if canon:
            alias_of = canon
            key = _norm(canon)

    pan = _panel()
    hit = pan[pan["_key"] == key]
    if not hit.empty:
        r = hit.iloc[0]
        prob = float(r["ile_prob"])
        # backfill the approved-drug-set rank when the same agent is in it
        app = _approved()
        arow = app[app["_key"] == key]
        rank = int(arow["rank"].iloc[0]) if not arow.empty else None
        return {
            "query": name, "matched": r["name"], "alias_of": alias_of, "source": "panel",
            "class": _clean(r.get("cls")), "logD": _num(r, "logD"), "Vd": _num(r, "Vd"),
            "MW": _num(r, "MW"), "PPB": _num(r, "PPB"), "theory_score": _num(r, "theory_score"),
            "ile_prob": round(prob, 1), "category": _category(prob),
            "ad_flag": _clean(r.get("ad_flag")) or _ad_from_mw(_num(r, "MW")),
            "rank": rank, "n_ranked": len(app) if rank else None, "note": _clean(r.get("note")),
        }

    app = _approved()
    hit = app[app["_key"] == key]
    if not hit.empty:
        r = hit.iloc[0]
        prob = float(r["ile_prob"])
        return {
            "query": name, "matched": r["name"], "alias_of": alias_of, "source": "approved_drugs",
            "class": None, "logD": _num(r, "logD"), "Vd": _num(r, "Vd"),
            "MW": _num(r, "MW"), "PPB": _num(r, "PPB"), "theory_score": _num(r, "theory_score"),
            "ile_prob": round(prob, 1), "category": _clean(r.get("category")) or _category(prob),
            "ad_flag": _clean(r.get("ad_flag")) or _ad_from_mw(_num(r, "MW")),
            "rank": int(r["rank"]) if "rank" in r and not pd.isna(r["rank"]) else None,
            "n_ranked": len(app), "note": None,
        }
    return None


def suggest(name: str, k: int = 6) -> list[str]:
    """Best-effort suggestions when a name is not found (normalised substring match)."""
    key = _norm(name)
    if not key:
        return []
    token = key.split()[0]
    out, seen = [], set()
    for d in (_panel(), _approved()):
        m = d[d["_key"].str.contains(re.escape(token), na=False)]
        for n in m["name"].head(k):
            if n.lower() not in seen:
                seen.add(n.lower())
                out.append(n)
    return out[:k]


def count_bundled() -> int:
    """Number of distinct bundled agents available for offline lookup."""
    keys = set(_panel()["_key"]) | set(_approved()["_key"])
    return len(keys)
