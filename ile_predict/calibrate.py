"""Data-driven calibration: logistic model of ILE responsiveness on logD.

Design note. Volume of distribution is intentionally EXCLUDED from the primary
model. In the real drug space, Vd is strongly collinear with lipophilicity
(lipophilic bases already have high Vd), so a naive fit flips its coefficient
positive -- the opposite of what the lipid-sink theory predicts. logD alone is
the robust, theory-consistent discriminator (LOO-CV AUC ~0.93).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

_DATA = Path(__file__).resolve().parent.parent / "data"


class CalibratedILEModel:
    """Logistic model P(responder | logD), fit on the reference label set."""

    def __init__(self):
        self._mean = None
        self._std = None
        self._coef = None
        self._intercept = None

    def fit(self, logd: np.ndarray, y: np.ndarray) -> "CalibratedILEModel":
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        logd = np.asarray(logd, float).reshape(-1, 1)
        scaler = StandardScaler().fit(logd)
        clf = LogisticRegression(C=0.5, max_iter=1000).fit(scaler.transform(logd), y)
        self._mean, self._std = scaler.mean_[0], scaler.scale_[0]
        self._coef, self._intercept = clf.coef_[0][0], clf.intercept_[0]
        return self

    def predict_proba(self, logd) -> np.ndarray:
        z = self._coef * ((np.asarray(logd, float) - self._mean) / self._std) + self._intercept
        return 1.0 / (1.0 + np.exp(-z))

    @classmethod
    def from_reference(cls, panel_csv: str | Path = None,
                       labels_csv: str | Path = None) -> "CalibratedILEModel":
        """Fit using the bundled reference labels joined to a scored panel."""
        panel = pd.read_csv(panel_csv or _DATA / "panel_scored.csv")
        labels = pd.read_csv(labels_csv or _DATA / "reference_labels.csv")
        logd_col = "logD" if "logD" in panel.columns else "Lipophilicity_AstraZeneca"
        m = panel.merge(labels, on="name", how="inner").dropna(subset=[logd_col, "ile_label"])
        return cls().fit(m[logd_col].values, m["ile_label"].astype(int).values)

    def loo_cv_auc(self, logd: np.ndarray, y: np.ndarray) -> float:
        """Leave-one-out cross-validated ROC-AUC (honest small-sample estimate)."""
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import LeaveOneOut
        from sklearn.metrics import roc_auc_score
        X = np.asarray(logd, float).reshape(-1, 1)
        y = np.asarray(y, int)
        p = np.zeros(len(y))
        for tr, te in LeaveOneOut().split(X):
            sc = StandardScaler().fit(X[tr])
            clf = LogisticRegression(C=0.5, max_iter=1000).fit(sc.transform(X[tr]), y[tr])
            p[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
        return float(roc_auc_score(y, p))
