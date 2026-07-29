"""Data-driven calibration: logistic model of ILE responsiveness on logD.

Design note. Volume of distribution is intentionally EXCLUDED from the primary
model. In the real drug space, Vd is strongly collinear with lipophilicity
(lipophilic bases already have high Vd), so a naive fit flips its coefficient
positive -- the opposite of what the lipid-sink theory predicts. logD alone is
the robust, theory-consistent discriminator (LOO-CV AUC ~0.93).
"""
from __future__ import annotations
import json
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
        self._train_logd = None
        self._train_y = None

    def fit(self, logd: np.ndarray, y: np.ndarray) -> "CalibratedILEModel":
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        logd = np.asarray(logd, float).reshape(-1, 1)
        y = np.asarray(y, int)
        scaler = StandardScaler().fit(logd)
        clf = LogisticRegression(C=0.5, max_iter=1000).fit(scaler.transform(logd), y)
        self._mean, self._std = scaler.mean_[0], scaler.scale_[0]
        self._coef, self._intercept = clf.coef_[0][0], clf.intercept_[0]
        self._train_logd = logd.ravel()
        self._train_y = y
        return self

    def predict_proba(self, logd) -> np.ndarray:
        z = self._coef * ((np.asarray(logd, float) - self._mean) / self._std) + self._intercept
        return 1.0 / (1.0 + np.exp(-z))

    def predict_proba_interval(
        self,
        logd,
        *,
        alpha: float = 0.05,
        n_bootstrap: int = 300,
        random_state: int = 42,
    ) -> dict:
        """Bootstrap confidence interval for probability predictions."""
        if self._train_logd is None or self._train_y is None:
            raise ValueError("Model must be fit before requesting intervals.")
        x = np.asarray(logd, float)
        point = self.predict_proba(x)
        rng = np.random.default_rng(random_state)
        boots = []
        n = len(self._train_y)
        attempts = 0
        max_attempts = n_bootstrap * 5
        while len(boots) < n_bootstrap and attempts < max_attempts:
            attempts += 1
            idx = rng.integers(0, n, size=n)
            ys = self._train_y[idx]
            if len(np.unique(ys)) < 2:
                continue
            m = CalibratedILEModel().fit(self._train_logd[idx], ys)
            boots.append(m.predict_proba(x))
        if not boots:
            lo = point
            hi = point
            used = 0
        else:
            boot_arr = np.asarray(boots, dtype=float)
            lo = np.quantile(boot_arr, alpha / 2.0, axis=0)
            hi = np.quantile(boot_arr, 1.0 - alpha / 2.0, axis=0)
            used = len(boots)
        return {
            "point": point,
            "lower": lo,
            "upper": hi,
            "level": 1.0 - alpha,
            "n_bootstrap": used,
        }

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

    @staticmethod
    def stability_auc(logd: np.ndarray, y: np.ndarray, random_state: int = 42) -> dict:
        """Calibration stability from LOO and repeated stratified CV."""
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import roc_auc_score
        from sklearn.model_selection import LeaveOneOut, RepeatedStratifiedKFold, StratifiedKFold
        from sklearn.preprocessing import StandardScaler

        X = np.asarray(logd, float).reshape(-1, 1)
        y = np.asarray(y, int)

        loo_p = np.zeros(len(y))
        for tr, te in LeaveOneOut().split(X):
            sc = StandardScaler().fit(X[tr])
            clf = LogisticRegression(C=0.5, max_iter=1000).fit(sc.transform(X[tr]), y[tr])
            loo_p[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
        loo_auc = float(roc_auc_score(y, loo_p))

        def _fold_scores(splitter):
            scores = []
            for tr, te in splitter.split(X, y):
                sc = StandardScaler().fit(X[tr])
                clf = LogisticRegression(C=0.5, max_iter=1000).fit(sc.transform(X[tr]), y[tr])
                p = clf.predict_proba(sc.transform(X[te]))[:, 1]
                scores.append(float(roc_auc_score(y[te], p)))
            return scores

        k = max(3, min(5, int(np.bincount(y).min())))
        skf_scores = _fold_scores(
            StratifiedKFold(n_splits=k, shuffle=True, random_state=random_state)
        )
        rskf_scores = _fold_scores(
            RepeatedStratifiedKFold(
                n_splits=k, n_repeats=10, random_state=random_state
            )
        )
        return {
            "loo_auc": loo_auc,
            "stratified_kfold_auc_mean": float(np.mean(skf_scores)),
            "stratified_kfold_auc_std": float(np.std(skf_scores)),
            "repeated_kfold_auc_mean": float(np.mean(rskf_scores)),
            "repeated_kfold_auc_std": float(np.std(rskf_scores)),
            "n_splits": int(k),
            "n_repeats": 10,
        }


def reference_label_metadata(path: str | Path = None) -> dict:
    """Read metadata for the bundled reference label set."""
    p = Path(path or (_DATA / "reference_labels.metadata.json"))
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as f:
        return json.load(f)


def reference_performance_report(
    panel_csv: str | Path = None,
    labels_csv: str | Path = None,
    min_group_n: int = 8,
) -> dict:
    """Performance report with overall and subgroup metrics."""
    from sklearn.metrics import roc_auc_score

    panel = pd.read_csv(panel_csv or _DATA / "panel_scored.csv")
    labels = pd.read_csv(labels_csv or _DATA / "reference_labels.csv")
    logd_col = "logD" if "logD" in panel.columns else "Lipophilicity_AstraZeneca"
    cls_col = "class" if "class" in panel.columns else ("cls" if "cls" in panel.columns else None)
    m = panel.merge(labels, on="name", how="inner").dropna(subset=[logd_col, "ile_label"])
    y = m["ile_label"].astype(int).values
    logd = m[logd_col].astype(float).values
    model = CalibratedILEModel().fit(logd, y)
    pred = model.predict_proba(logd)
    overall_auc = float(roc_auc_score(y, pred))
    m["_mw_band"] = pd.cut(
        m["MW"].astype(float),
        bins=[-np.inf, 100, 200, 400, 600, np.inf],
        labels=["<100", "100-200", "200-400", "400-600", ">600"],
    )
    ad_series = m["AD"] if "AD" in m.columns else (m["ad_flag"] if "ad_flag" in m.columns else pd.Series(["unknown"] * len(m)))
    m["_ad_flag_norm"] = ad_series.fillna("unknown")

    def _group_auc(frame: pd.DataFrame, by_col: str):
        out = []
        for key, g in frame.groupby(by_col):
            yy = g["ile_label"].astype(int).values
            if len(g) < min_group_n or len(np.unique(yy)) < 2:
                continue
            pp = model.predict_proba(g[logd_col].astype(float).values)
            out.append(
                {
                    "group": str(key),
                    "n": int(len(g)),
                    "responders": int((yy == 1).sum()),
                    "non_responders": int((yy == 0).sum()),
                    "auc": float(roc_auc_score(yy, pp)),
                }
            )
        return sorted(out, key=lambda x: (-x["n"], x["group"]))

    return {
        "metadata": reference_label_metadata(),
        "n_reference": int(len(m)),
        "overall_auc_in_sample": overall_auc,
        "stability": CalibratedILEModel.stability_auc(logd, y),
        "subgroups": {
            "class": _group_auc(m, cls_col) if cls_col else [],
            "mw_band": _group_auc(m, "_mw_band"),
            "applicability": _group_auc(m, "_ad_flag_norm"),
        },
    }
