"""
ile-predict: Lipid emulsion (ILE) amenability prediction from chemical structure.

A research prototype that predicts how amenable a toxic substance is to
intravenous lipid emulsion ("lipid sink") therapy, using ADMET-AI physicochemical
predictions (logD, volume of distribution) and a logistic model calibrated against
a literature-consensus ILE-evidence reference set.

NOT A CLINICAL DECISION TOOL. Hypothesis generation for research only.
"""
__version__ = "0.2.0"

# Light-weight, always-available API (offline lookup + calibration).
from .lookup import lookup, suggest, count_bundled       # noqa: F401
from .calibrate import (                                 # noqa: F401
    CalibratedILEModel,
    reference_label_metadata,
    reference_performance_report,
)

# Structure scoring needs the optional [full] stack (ADMET-AI, RDKit). The imports
# below are lazy inside the functions, so importing this package stays light; a
# missing stack only raises when you actually score a novel structure.
from .predict import score_smiles, score_batch           # noqa: F401
