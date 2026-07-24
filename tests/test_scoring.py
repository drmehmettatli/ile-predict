"""Unit tests for the theory score and calibration model (no network / no ADMET-AI)."""
import numpy as np
from ile_predict.predict import theory_score, category
from ile_predict.calibrate import CalibratedILEModel


def test_theory_score_monotonic_in_logd():
    # higher effective lipophilicity -> higher score (Vd fixed)
    assert theory_score(4.0, 2.0) > theory_score(1.0, 2.0)


def test_theory_score_penalizes_high_vd():
    # same logD, higher Vd -> lower score
    assert theory_score(3.0, 1.0) > theory_score(3.0, 30.0)


def test_theory_score_bounds():
    for d in (-2, 0, 3, 8):
        for v in (0.1, 1, 50):
            s = theory_score(d, v)
            assert 0.0 <= s <= 100.0


def test_category_thresholds():
    assert category(80) == "High"
    assert category(50) == "Medium"
    assert category(10) == "Low"


def test_calibration_separates_classes():
    # synthetic: responders have high logD, non-responders low
    logd = np.array([3.5, 3.0, 2.8, 0.2, -0.5, 0.0])
    y = np.array([1, 1, 1, 0, 0, 0])
    m = CalibratedILEModel().fit(logd, y)
    p_hi = m.predict_proba([3.5])[0]
    p_lo = m.predict_proba([-0.5])[0]
    assert p_hi > 0.5 > p_lo


def test_loo_cv_auc_perfect_on_separable():
    logd = np.array([3.5, 3.2, 3.0, 2.8, 0.3, 0.1, -0.2, -0.5])
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0])
    auc = CalibratedILEModel().loo_cv_auc(logd, y)
    assert auc >= 0.9
