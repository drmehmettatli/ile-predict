"""Command-line interface for ile-predict.

Offline-first: named drugs are answered instantly from the bundled tables with no
heavy dependencies. A raw SMILES, or a name not in the tables, is scored with the
full ADMET-AI stack when it is installed (``pip install 'ile-predict[full]'``).

Examples
--------
    ile-predict "lurasidone, nimodipine, oleandrin, metformin"
    ile-predict "verapamil" --json
    ile-predict --smiles "CCN(CC)CC(=O)Nc1c(C)cccc1C" --out results.csv
"""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .lookup import lookup, suggest, count_bundled

_COLS = ["query", "matched", "source", "class", "logD", "Vd", "MW",
         "theory_score", "ile_prob", "ci_low", "ci_high", "category", "ad_flag",
         "confidence_policy", "status", "rank", "note"]

_DISCLAIMER = "[!] Research/education prototype. NOT a clinical decision tool."
_COLAB = "https://colab.research.google.com/github/drmehmettatli/ile-predict/blob/main/notebooks/ile_predict_colab.ipynb"


def _admet_available() -> bool:
    import importlib.util
    return all(importlib.util.find_spec(m) for m in ("admet_ai", "rdkit"))


def _score_novel(pairs):
    """pairs: list of (query_label, smiles). Returns list of record dicts.

    Imports the heavy stack lazily so the offline path never pays for it.
    """
    from .data import canonicalize
    from .predict import score_batch, category
    from .calibrate import CalibratedILEModel

    labels = [p[0] for p in pairs]
    smis = [canonicalize(p[1]) for p in pairs]
    keep = [(lab, s) for lab, s in zip(labels, smis) if s]
    if not keep:
        return []
    df = score_batch([s for _, s in keep])
    model = CalibratedILEModel.from_reference()
    df["ile_prob"] = (model.predict_proba(df["logD"].values) * 100).round(1)
    recs = []
    for (lab, _), (_, row) in zip(keep, df.iterrows()):
        prob = float(row["ile_prob"])
        ad_flag = row["ad_flag"]
        recs.append({
            "query": lab, "matched": lab, "source": "computed",
            "class": None, "logD": round(float(row["logD"]), 3),
            "Vd": round(float(row["Vd"]), 3), "MW": round(float(row["MW"]), 1),
            "theory_score": round(float(row["theory_score"]), 1), "ile_prob": round(prob, 1),
            "category": category(prob), "ad_flag": ad_flag, "rank": None, "note": None,
            "status": "computed-from-structure",
            "confidence_policy": "low" if "domain-edge" in str(ad_flag).lower() else "standard",
        })
    return recs


def _resolve_names_to_smiles(names):
    from .data import name_to_smiles
    return [(n, name_to_smiles(n)) for n in names]


def _print_table(rows):
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        for c in _COLS:
            if c not in df.columns:
                df[c] = None
        with pd.option_context("display.width", 200, "display.max_colwidth", 40):
            print(df[_COLS].to_string(index=False))
    except Exception:
        for r in rows:
            print(r)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="ile-predict",
        description="Predict lipid emulsion (ILE) amenability from structure. "
                    "RESEARCH ONLY, not clinical advice.")
    ap.add_argument("query", nargs="?", help="comma-separated drug names")
    ap.add_argument("--smiles", help="comma-separated SMILES (forces full-stack scoring)")
    ap.add_argument("--out", help="write full results to CSV")
    ap.add_argument("--json", action="store_true", help="print results as JSON")
    ap.add_argument("--no-ci", action="store_true",
                    help="skip bootstrap confidence intervals around calibrated probability")
    ap.add_argument("--no-compute", action="store_true",
                    help="offline only: never fall back to ADMET-AI for unknown names")
    ap.add_argument("--version", action="version", version=f"ile-predict {__version__}")
    args = ap.parse_args(argv)

    rows, unresolved = [], []

    if args.smiles:
        smis = [s.strip() for s in args.smiles.split(",") if s.strip()]
        if not _admet_available():
            print(_novel_help(smis, kind="SMILES"), file=sys.stderr)
            return 2
        rows = _score_novel([(s, s) for s in smis])

    elif args.query:
        names = [n.strip() for n in args.query.split(",") if n.strip()]
        missing = []
        for n in names:
            rec = lookup(n)
            if rec:
                rows.append(rec)
            else:
                missing.append(n)
        if missing and not args.no_compute and _admet_available():
            resolved = _resolve_names_to_smiles(missing)
            ok = [(n, s) for n, s in resolved if s]
            unresolved = [n for n, s in resolved if not s]
            rows += _score_novel(ok)
        elif missing:
            unresolved = missing
    else:
        ap.error("provide drug names or --smiles")

    if not rows and unresolved:
        print(_novel_help(unresolved, kind="name"), file=sys.stderr)
        return 1

    _attach_confidence(rows, enabled=not args.no_ci)

    if args.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        _print_table(rows)
        print("\n" + _DISCLAIMER)

    for n in unresolved:
        sug = suggest(n)
        tip = f"  did you mean: {', '.join(sug)}" if sug else ""
        print(f"[not found] {n}{tip}", file=sys.stderr)
    _print_reliability_notes(rows)

    if args.out:
        try:
            import pandas as pd
            pd.DataFrame(rows).to_csv(args.out, index=False)
            print(f"written: {args.out}")
        except Exception as e:
            print(f"could not write {args.out}: {e}", file=sys.stderr)
    return 0


def _attach_confidence(rows, enabled=True):
    if not enabled or not rows:
        return
    try:
        from .calibrate import CalibratedILEModel
        model = CalibratedILEModel.from_reference()
        to_score = [r for r in rows if r.get("logD") is not None]
        if not to_score:
            return
        interval = model.predict_proba_interval([r["logD"] for r in to_score], n_bootstrap=150)
        for i, r in enumerate(to_score):
            r["ci_low"] = round(float(interval["lower"][i] * 100.0), 1)
            r["ci_high"] = round(float(interval["upper"][i] * 100.0), 1)
    except Exception:
        return


def _novel_help(items, kind):
    n = count_bundled()
    head = f"None of the requested {kind}(s) are in the {n} bundled agents, "
    return (head + "and the full prediction stack is not installed.\n"
            "  To score novel structures locally:  pip install 'ile-predict[full]'\n"
            f"  Or run with zero local setup in Colab:\n    {_COLAB}")


def _print_reliability_notes(rows):
    for r in rows:
        if r.get("confidence_policy") == "low":
            label = r.get("matched") or r.get("query") or "item"
            print(
                f"[low confidence] {label}: outside ADMET-AI domain; treat probability as uncertain.",
                file=sys.stderr,
            )


if __name__ == "__main__":
    raise SystemExit(main())
