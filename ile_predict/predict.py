"""ADMET-AI prediction + mechanistic (theory-based) lipid-sink score.

The theory score encodes the lipid-sink hypothesis:
  * lipophilicity is a NECESSARY condition (gate) -> uses logD7.4 (effective,
    ionization-aware lipophilicity), not raw logP;
  * volume of distribution MODULATES it (low Vd = drug stays central = easier
    to sequester).
"""
from __future__ import annotations
from typing import Iterable
import numpy as np
import pandas as pd

from .data import applicability_flag

_MODEL = None


def _model():
    global _MODEL
    if _MODEL is None:
        from admet_ai import ADMETModel
        _MODEL = ADMETModel()
    return _MODEL


def theory_score(logd: float, vd: float) -> float:
    """Mechanistic lipid-sink amenability score in [0, 100]."""
    vd = max(vd, 0.2)
    lipo = np.clip(logd / 3.0, 0, 1) * 100.0            # logD 0->0, 3->100
    vd_mod = np.clip(1 - np.log10(vd), 0, 1)            # Vd 1->1.0, 10->0
    return float(round(lipo * (0.6 + 0.4 * vd_mod), 1))


def category(prob: float) -> str:
    return "High" if prob >= 70 else ("Medium" if prob >= 40 else "Low")


def score_batch(smiles: Iterable[str]) -> pd.DataFrame:
    """Run ADMET-AI on a list of SMILES and attach the theory score.

    Returns a DataFrame with logP, logD, Vd, PPB, MW, theory_score, ad_flag.
    """
    smiles = list(smiles)
    preds = _model().predict(smiles=smiles).rename_axis("smiles").reset_index()
    preds["logD"] = preds["Lipophilicity_AstraZeneca"]
    preds["Vd"] = preds["VDss_Lombardo"]
    preds["PPB"] = preds["PPBR_AZ"].clip(0, 100)
    preds["MW"] = preds["molecular_weight"]
    preds["theory_score"] = [
        theory_score(d, v) for d, v in zip(preds["logD"], preds["Vd"])
    ]
    preds["ad_flag"] = preds["MW"].map(applicability_flag)
    cols = ["smiles", "logP", "logD", "Vd", "PPB", "MW", "theory_score", "ad_flag"]
    return preds[cols]


def score_smiles(smiles: str) -> dict:
    """Convenience wrapper for a single SMILES."""
    return score_batch([smiles]).iloc[0].to_dict()
