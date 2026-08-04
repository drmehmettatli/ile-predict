# Methods

## 1. Inputs
Compound names are resolved to isomeric SMILES via **PubChem** (`pubchempy`) and
canonicalised with **RDKit**. Physicochemical / ADMET descriptors are predicted by
**ADMET-AI** (Chemprop graph neural networks trained on Therapeutics Data Commons):
`logP`, `logD7.4` (`Lipophilicity_AstraZeneca`), steady-state volume of distribution
`Vd` (`VDss_Lombardo`, L/kg), and plasma protein binding (`PPBR_AZ`).

## 2. Mechanistic (theory) score
The lipid-sink hypothesis: an intravascular lipid phase sequesters drug in
proportion to its *effective* lipophilicity, and low-Vd drugs (which stay central)
are easier to recapture.

```
Lipo   = clip(logD / 3, 0, 1) * 100            # necessary condition (gate)
VdMod  = clip(1 - log10(max(Vd, 0.2)), 0, 1)   # Vd 1 L/kg -> 1.0 ; >=10 -> 0
score  = Lipo * (0.6 + 0.4 * VdMod)            # Vd modulates 60-100% of Lipo
```
`logD7.4` is used instead of raw `logP` because ionisation governs the neutral
fraction available to partition; acidic/highly-bound drugs (e.g. valsartan
logP 4.2 but logD -0.03) are otherwise mis-ranked.

## 3. Data-driven calibration
A logistic regression is fit on a reference set of drugs labelled as ILE
*responders* / *non-responders* from literature consensus (`data/reference_labels.csv`).
Performance is reported as **leave-one-out cross-validated ROC-AUC** (honest for the
small sample): logD-only ≈ 0.93; raw logP ≈ 0.90; the mechanistic theory score ≈ 0.92.

### Why volume of distribution is excluded from the primary model
Adding `log10(Vd)` raises AUC slightly (~0.95) but its fitted coefficient is
**positive**, opposite to theory, which expects *low* Vd to favour ILE. This is a
collinearity artefact: in the real drug space, lipophilic bases already have high
Vd, so Vd tracks lipophilicity rather than acting as an independent recapture term.
The primary calibrated model therefore uses **logD alone**; the theory score retains
the (mechanistically correct) Vd penalty and is reported alongside as a prior.

### What a one-variable calibrated model does and does not add
Because the calibrated model has a single predictor, the logistic link is a monotone
transformation of logD. The **ranking** it produces is therefore identical to ranking
compounds by logD; calibration adds an interpretable probability scale and a
data-anchored threshold, not new discriminative information. Any claim that the
calibrated score orders compounds better than logD alone would be false. The only
component that re-orders compounds relative to logD is the mechanistic theory score,
through its Vd modulator.

## 4. Applicability domain
ADMET-AI is trained on drug-like molecules. Predictions for MW <100 (gases, small
solvents) or MW >600 (large natural products, depot esters) are flagged
`domain-edge` and treated as unreliable; a low score there means "uncertain", not
"not a candidate" (e.g. aconitine, veratridine).

## 5. Mechanism caveats (score ≠ clinical benefit)
High lipophilicity does not imply ILE benefit when the toxic mechanism is not
reversible by lowering free parent drug: carbon monoxide (carboxyhaemoglobin),
vitamin D (downstream hypercalcaemia via calcitriol), aliphatic hydrocarbons
(aspiration pneumonitis, a local injury). Volatile hydrocarbons additionally have
distinct redistribution/exhalation kinetics.

## 6. Known limitations

**Reference labels.** Labels are expert-assigned literature consensus, not a single
validated gold standard; residual false positives (e.g. warfarin, digoxin, high logD
but specific antidotes) are handled by a clinical-note layer, not the score.
Calibration to measured *in-vitro* sequestration data is the recommended next step.

**Negative volume-of-distribution predictions.** `VDss_Lombardo` is an unbounded
regression head and returns physiologically impossible negative values for 48 of the
209 panel agents, that is 23.0 per cent. All 48 lie inside the molecular-weight
applicability domain, so this is not an out-of-domain effect. The failures are
systematic rather than random: affected agents have lower predicted logD (median 0.63
versus 2.00, Mann-Whitney p = 3.1e-08), lower molecular weight (median 226.8 versus
295.4 Da) and lower predicted plasma protein binding (median 52.8 versus 69.5 per
cent). The most extreme cases are oxcarbazepine (-5.90 L/kg), phenytoin (-4.78),
carbamazepine (-4.53), gabapentin (-4.53), glyburide (-4.44) and naproxen (-4.14),
all of which are genuinely low-Vd drugs. The model therefore appears to rank them
correctly and then overshoot the physiological floor; plasma volume alone is roughly
0.04 L/kg. The artefact looks monotone, so rank order may survive while the scale
does not. The theory score guards against it with `max(Vd, 0.2)`, which is
rank-preserving but crude, and the affected compounds should be treated as having no
usable Vd estimate rather than a low one.
