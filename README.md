# ile-predict

**The TATLI framework — Toxin Amenability To Lipid Infusion.**
**Predicting intravenous lipid emulsion (ILE / "lipid rescue") amenability from chemical structure.**

Developed in the **Van Computational Emergency Toxicology** group at Van Yüzüncü Yıl
University and the University of Health Sciences (Van Training and Research Hospital),
Van, Türkiye, by **Mehmet Tatlı** ([ORCID 0000-0001-5907-9161](https://orcid.org/0000-0001-5907-9161)).

`ile-predict` estimates how amenable a toxic substance is to intravenous lipid
emulsion therapy, directly from its structure. It combines
[ADMET-AI](https://github.com/swansonk14/admet_ai) physicochemical predictions
(effective lipophilicity `logD7.4`, volume of distribution) with a logistic model
**calibrated against a literature-consensus ILE-evidence reference set**.

> ⚠️ **Research / education prototype — NOT a clinical decision tool.**
> The score reflects only *physicochemical* amenability under the lipid-sink
> hypothesis. It does **not** account for clinical severity, available specific
> antidotes, or whether the toxin's mechanism is even reversible by lowering free
> drug (e.g. carbon-monoxide/haemoglobin binding, vitamin-D hypercalcaemia). Real
> ILE decisions follow current guidelines (e.g. ASRA/AAGBI) and clinical judgment.

## Which path is for you?

| You want to… | Use | Install |
|---|---|---|
| **Look up a named drug** (score, category, rank, caveats) | the **search dashboard** — `app/dashboard.html`, open in any browser | none |
| Look up / batch drugs on the command line or in Python | the **light** package (`ile-predict "lurasidone, nimodipine"`) | `pip install "ile-predict @ git+https://github.com/drmehmettatli/ile-predict.git"` |
| **Score a novel structure** (a SMILES or a drug not yet bundled) | the **full** stack, or the notebook below | `pip install "ile-predict[full] @ git+https://github.com/drmehmettatli/ile-predict.git"` |
| Score novel structures with **zero local install** | the Colab notebook → | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/drmehmettatli/ile-predict/blob/main/notebooks/ile_predict_colab.ipynb) |

The dashboard and the light package answer instantly for **~2,900 pre-scored agents**
(209-agent curated panel + 2,845 approved drugs) **with no heavy dependencies** — the
common "what does *this drug* score?" question needs no ADMET-AI, no PyTorch, no
download. Only a genuinely novel structure invokes the full prediction stack.

## What it does

1. **Structure in** — a drug/toxin name (resolved to SMILES via PubChem) or a raw SMILES.
2. **Predict** — ADMET-AI computes `logP`, `logD7.4`, `Vd`, plasma protein binding.
3. **Score** — a mechanistic *theory score* (lipophilicity gate × Vd modulator) **and**
   a *calibrated probability* from the logistic model, reported as the **TATLI amenability
   score** (Toxin Amenability To Lipid Infusion).
4. **Flag** — molecules outside ADMET-AI's drug-like domain (MW <100 or >600) are
   marked `domain-edge`; a low score there means "unreliable", not "not a candidate".

## Quick start

**Zero install — the dashboard.** Download `app/dashboard.html` and open it in a
browser. Type a drug name in the search box to get its amenability score, category,
rank among 2,845 approved drugs, applicability-domain flag, and mechanism caveats.
Common abbreviations and brand/street names resolve automatically (e.g. `THC` →
dronabinol, `CBD` → cannabidiol, `ASA` → aspirin, `Seroquel` → quetiapine); the full
map is `data/aliases.csv`.

**Light install — instant offline lookup** (core only: pandas, numpy, scikit-learn):

```bash
pip install "ile-predict @ git+https://github.com/drmehmettatli/ile-predict.git"
# or from a clone:  pip install -r requirements.txt && pip install -e .

ile-predict "lurasidone, nimodipine, oleandrin, metformin"   # answered from bundled tables
ile-predict "verapamil" --json                               # machine-readable output
```

```python
from ile_predict import lookup

rec = lookup("lurasidone")
print(rec["ile_prob"], rec["category"], rec["rank"])   # 99.2 High 66
```

**Full install — score novel structures** (adds ADMET-AI + RDKit; first run downloads
the ADMET-AI models):

```bash
pip install "ile-predict[full] @ git+https://github.com/drmehmettatli/ile-predict.git"
# or from a clone:  pip install -r requirements-full.txt && pip install -e '.[full]'

# a SMILES, or any drug not in the bundled tables, is computed on the fly:
ile-predict --smiles "CCN(CC)CC(=O)Nc1c(C)cccc1C" --out results.csv
```

```python
from ile_predict import score_smiles
from ile_predict.calibrate import CalibratedILEModel

rec = score_smiles("CCN(CC)CC(=O)Nc1c(C)cccc1C")   # lidocaine
model = CalibratedILEModel.from_reference()
prob = model.predict_proba([rec["logD"]])[0] * 100
```

## Method (short)

The lipid-sink hypothesis says a circulating lipid phase sequesters drug in
proportion to its effective lipophilicity. We therefore gate on `logD7.4` (not raw
`logP`, which mis-ranks acidic/highly-bound drugs), and modulate by volume of
distribution. The **calibrated model uses logD alone**: adding Vd nudges AUC up but
flips its coefficient positive — an artefact of Vd/lipophilicity collinearity in
real drugs — so it is excluded from the primary model. See [`docs/methods.md`](docs/methods.md).

- Reference label set: `data/reference_labels.csv` (76 drugs, expert-assigned from
  ILE literature consensus — **replace/extend with a systematic review for production use**).
- Scored compound panel: `data/panel_scored.csv` (209 agents across 31 classes).
- Approved-drug screen: `data/approved_drugs_scored.csv` (2,845 drugs, scored and
  ranked — the set used in the manuscript; powers offline lookup).
- Interactive search dashboard: `app/dashboard.html` (self-contained, open in a browser).
- Reported performance: leave-one-out cross-validated ROC-AUC ≈ **0.93**.

## Roadmap

- Expand the calibration label set via a systematic review of ILE case reports.
- Recalibrate against measured *in-vitro* lipid:aqueous sequestration data.
- Package a hosted web app (Streamlit/FastAPI) and a DrugBank/ChEMBL-scale screen.

## How to cite

If you use this software or the TATLI amenability score, please cite the archived
release and the accompanying manuscript:

> Tatlı M. *ile-predict (TATLI framework): a structure-based prediction system for
> intravenous lipid emulsion amenability in acute poisoning.* Van Computational
> Emergency Toxicology, Van Yüzüncü Yıl University and University of Health Sciences,
> Van, Türkiye. Zenodo. https://doi.org/10.5281/zenodo.21547460 (concept DOI; use the
> versioned DOI for the exact release you ran).

> Tatlı M. *Physicochemical amenability to intravenous lipid emulsion does not predict
> which poisonings have been studied: a structure-based screen of 209 toxic agents and
> 2845 approved drugs.* Manuscript in preparation, 2026.

ADMET predictions: Swanson et al., *ADMET-AI*, **Bioinformatics** 2024;40(7):btae416.
Structures: PubChem. Calibration labels: literature consensus (expert-assigned).
See `CITATION.cff`. Licensed under MIT (`LICENSE`).

## Disclaimer

This software is provided for research and educational purposes only and does not
constitute medical advice. No warranty is made as to accuracy or fitness for
clinical use. Do not use it to guide patient care.
