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
**positive** — opposite to theory, which expects *low* Vd to favour ILE. This is a
collinearity artefact: in the real drug space, lipophilic bases already have high
Vd, so Vd tracks lipophilicity rather than acting as an independent recapture term.
The primary calibrated model therefore uses **logD alone**; the theory score retains
the (mechanistically correct) Vd penalty and is reported alongside as a prior.

## 4. Applicability domain
ADMET-AI is trained on drug-like molecules. Predictions for MW <100 (gases, small
solvents) or MW >600 (large natural products, depot esters) are flagged
`domain-edge` and treated as unreliable — a low score there means "uncertain", not
"not a candidate" (e.g. aconitine, veratridine).

## 5. Mechanism caveats (score ≠ clinical benefit)
High lipophilicity does not imply ILE benefit when the toxic mechanism is not
reversible by lowering free parent drug: carbon monoxide (carboxyhaemoglobin),
vitamin D (downstream hypercalcaemia via calcitriol), aliphatic hydrocarbons
(aspiration pneumonitis, a local injury). Volatile hydrocarbons additionally have
distinct redistribution/exhalation kinetics.

## 6. Known limitations
Labels are expert-assigned literature consensus, not a single validated gold
standard; residual false positives (e.g. warfarin, digoxin — high logD but specific
antidotes) are handled by a clinical-note layer, not the score. Calibration to
measured *in-vitro* sequestration data is the recommended next step.
