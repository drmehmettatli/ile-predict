"""Reproducible enumeration of the two drug classes used for the within-class
comparisons reported in the manuscript (Results, "Predicted amenability is
discordant with published evidence"; Table 2; Figure 4).

Dihydropyridine calcium channel blockers are enumerated structurally, by an
RDKit substructure match to the 1,4-dihydropyridine ring, so that the class
membership does not depend on a hand-written name list. Two post-filters are
applied and both are declared in the manuscript: the endogenous cofactor NADH
matches the ring but is not a calcium channel blocker and is dropped, and
levamlodipine is the resolved S-enantiomer of amlodipine and is counted with
it. The result is 13 distinct agents.

Antipsychotics cannot be captured by a single substructure, because the class
spans phenothiazines, butyrophenones, thioxanthenes, benzisoxazoles and
several other scaffolds. They are therefore enumerated by generic name against
the established typical and atypical agents. The result is 32 agents.

Usage:
    python -m ile_predict.class_enumeration path/to/drugbank_ile_scored.csv
"""

from __future__ import annotations

import sys

import pandas as pd

DHP_SMARTS = "[#6]1=[#6][#6][#6]=[#6][#7]1"

DHP_EXCLUDE = {"nadh"}
DHP_MERGE = {"levamlodipine": "amlodipine"}

ANTIPSYCHOTICS = [
    "chlorpromazine", "fluphenazine", "perphenazine", "prochlorperazine",
    "trifluoperazine", "thioridazine", "mesoridazine", "haloperidol",
    "droperidol", "pimozide", "loxapine", "molindone", "thiothixene",
    "flupentixol", "zuclopenthixol", "clozapine", "olanzapine", "quetiapine",
    "risperidone", "paliperidone", "ziprasidone", "aripiprazole",
    "brexpiprazole", "cariprazine", "asenapine", "iloperidone", "lurasidone",
    "amisulpride", "sulpiride", "sertindole", "pimavanserin", "lumateperone",
]


def dihydropyridines(df: pd.DataFrame) -> pd.DataFrame:
    """Return the dihydropyridine calcium channel blockers in ``df``.

    ``df`` must carry ``name``, ``smiles`` and ``ILE_prob`` columns.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    patt = Chem.MolFromSmarts(DHP_SMARTS)

    keep = []
    for _, row in df.iterrows():
        smiles = row.get("smiles")
        if not isinstance(smiles, str):
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or not mol.HasSubstructMatch(patt):
            continue
        name = str(row["name"]).lower()
        if name in DHP_EXCLUDE or name in DHP_MERGE:
            continue
        keep.append(row)

    out = pd.DataFrame(keep)
    return out.sort_values("ILE_prob", ascending=False).reset_index(drop=True)


def antipsychotics(df: pd.DataFrame) -> pd.DataFrame:
    """Return the antipsychotics in ``df``, ranked by calibrated probability."""
    names = df["name"].str.lower()
    out = df[names.isin(ANTIPSYCHOTICS)]
    return out.sort_values("ILE_prob", ascending=False).reset_index(drop=True)


def main(path: str) -> None:
    df = pd.read_csv(path)
    dhp = dihydropyridines(df)
    aps = antipsychotics(df)

    print(f"Dihydropyridine calcium channel blockers: n = {len(dhp)}")
    print(dhp[["name", "logD", "Vd_Lkg", "ILE_score", "ILE_prob"]].to_string(index=False))
    print()
    print(f"Antipsychotics: n = {len(aps)}")
    print(aps[["name", "logD", "Vd_Lkg", "ILE_score", "ILE_prob"]].to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
