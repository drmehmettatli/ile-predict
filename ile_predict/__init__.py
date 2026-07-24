"""
ile-predict: Lipid emulsion (ILE) amenability prediction from chemical structure.

A research prototype that predicts how amenable a toxic substance is to
intravenous lipid emulsion ("lipid sink") therapy, using ADMET-AI physicochemical
predictions (logD, volume of distribution) and a logistic model calibrated against
a literature-consensus ILE-evidence reference set.

NOT A CLINICAL DECISION TOOL. Hypothesis generation for research only.
"""
__version__ = "0.1.0"

from .predict import score_smiles, score_batch          # noqa: F401
from .calibrate import CalibratedILEModel               # noqa: F401
