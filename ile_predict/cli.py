"""Command-line interface:  ile-predict "verapamil,bupivacaine"  or  --smiles ... """
from __future__ import annotations
import argparse
import sys
import pandas as pd

from .data import name_to_smiles, canonicalize
from .predict import score_batch, category
from .calibrate import CalibratedILEModel


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Predict lipid emulsion (ILE) amenability. RESEARCH ONLY, not clinical advice.")
    ap.add_argument("query", nargs="?", help="comma-separated drug names")
    ap.add_argument("--smiles", help="comma-separated SMILES (skip name lookup)")
    ap.add_argument("--out", help="write results to CSV")
    args = ap.parse_args(argv)

    if args.smiles:
        names = args.smiles.split(",")
        smiles = [canonicalize(s.strip()) for s in names]
    elif args.query:
        names = [n.strip() for n in args.query.split(",")]
        smiles = [name_to_smiles(n) for n in names]
    else:
        ap.error("provide drug names or --smiles")

    valid = [(n, s) for n, s in zip(names, smiles) if s]
    if not valid:
        print("No resolvable structures.", file=sys.stderr)
        return 1

    df = score_batch([s for _, s in valid])
    df.insert(0, "query", [n for n, _ in valid])

    model = CalibratedILEModel.from_reference()
    df["ile_prob"] = (model.predict_proba(df["logD"].values) * 100).round(1)
    df["category"] = df["ile_prob"].map(category)

    show = df[["query", "logD", "Vd", "theory_score", "ile_prob", "category", "ad_flag"]]
    with pd.option_context("display.width", 160):
        print(show.round(2).to_string(index=False))
    print("\n[!] Research/education prototype. NOT a clinical decision tool.")
    if args.out:
        df.to_csv(args.out, index=False)
        print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
