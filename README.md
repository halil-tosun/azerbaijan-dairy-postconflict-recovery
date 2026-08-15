# Technical Efficiency, Productivity Growth, and Post-Conflict Recovery: The Fragility of Causal Claims in Azerbaijan's Dairy Sector, 2000-2024

## Reproducibility Package

This repository contains the complete reproducibility package
accompanying a manuscript that estimates technical efficiency,
decomposes productivity growth, and evaluates the effect of
Azerbaijan's post-2020 territorial reintegration on dairy-sector
technical efficiency, using a district-level panel (67 districts,
2000-2024). The empirical framework combines bootstrap Data Envelopment
Analysis (DEA), Stochastic Frontier Analysis (SFA), a bootstrap
Malmquist productivity index, a Simar-Wilson (2007) second-stage
determinants regression, and an interference-robust ("control-only
frontier") event-study design that addresses a Stable Unit Treatment
Value Assumption (SUTVA) concern in the simpler pooled-frontier
alternative.

This package is intentionally organized around the *study* rather than
any single journal submission. If the manuscript title, framing, or
target journal changes during peer review, this repository and its
contents remain valid without modification.

---

## Key Finding

Productivity growth in Azerbaijan's dairy sector over 2000-2024 was
driven predominantly by technological change rather than efficiency
convergence (positive in 11 of 13 economic regions; efficiency change
negative in all 13), and districts with larger dairy herds are
systematically less efficient despite mild increasing returns to
scale. For the causal question motivating the study -- did territorial
reintegration improve dairy-sector efficiency in reintegrated districts
-- a naive difference-in-differences estimate on pooled DEA scores
appears statistically significant (p = 0.002, classical SE); once
corrected for the treated districts' capacity to mechanically shift the
reference frontier against which every district is scored, and with
standard errors clustered at the district level, the estimated effect
is larger in magnitude (+0.116 vs. +0.057) but no longer statistically
distinguishable from zero (p = 0.300), and the design's own minimum
detectable effect (0.312, i.e., 74% of the sample mean efficiency
score) shows this is at least partly a question of statistical power
rather than a precisely estimated null. This is treated as a
substantive methodological finding in its own right, not a null result
to explain away, with direct relevance to a growing applied literature
that evaluates policy shocks using frontier-based efficiency measures
as outcome variables.

---

## Repository Overview

This repository follows open science and computational reproducibility
principles and includes:

- Complete Python source code for the full pipeline: data preparation;
  bootstrap DEA; stochastic frontier analysis; the second-stage
  determinants regression; bootstrap Malmquist productivity
  decomposition; the interference-robust control-frontier event-study
  design (clustered/robust standard errors, fractional-logit
  specification, minimum-detectable-effect and power calculations);
  augmented-technology sensitivity checks; robustness diagnostics; four
  additional robustness analyses added following external methodological
  review; and all main-text and Supporting Information figures
- Raw and analytical district-level datasets
- All main-text tables (Tables 1-9) and all Supporting Information
  tables (Tables S1-S18) reported in the manuscript
- All five main-text figures and two Supporting Information figures
- The manuscript's graphical abstract
- Comprehensive documentation of every model equation, data source, and
  the exact correspondence between each script and each reported
  statistic
- Software environment specifications

---

## Repository Structure

```text
Azerbaijan_Dairy_Replication/
├── code/
│   ├── 00_data_prep.py                       Data preparation and Table 1 sample-accounting check
│   ├── 01_dea_bootstrap.py                   Bootstrap DEA -- Table 2, Figure 2
│   ├── 02_sfa_models.py                      SFA models -- Table 3, SI Tables S1-S2
│   ├── 03_second_stage.py                    Second-stage determinants regression -- Table 5, SI S3.1/S5/S6/S9
│   ├── 04_malmquist_bootstrap.py             Bootstrap Malmquist -- Table 4, Figure 3, SI Figures S1-S2
│   ├── 05_event_study.py                     Legacy pooled-frontier event study (superseded by 09; retained for comparison)
│   ├── 06_labor_capital_land_robustness.py   Augmented-technology sensitivity checks -- Table 9, SI S10-S12
│   ├── 07_robustness_diagnostics.py          Balanced panel, VIF, district ranking -- SI S4/S7/S8
│   ├── 08_make_figures.py                    All main-text and SI figures (300 DPI) -- single authoritative figure script
│   ├── 08a_make_figure1_framework.py         Legacy/unused -- not referenced in the current manuscript; not called by run_all.py
│   ├── 09_control_frontier_event_study.py    Control-frontier event study, clustered/robust SE, fractional logit, MDE/power -- Tables 6-8, Figure 4
│   ├── 10_gpt_review_robustness.py           Additional robustness from external review -- SI S15-S18
│   └── run_all.py                            Runs the full pipeline in order
├── data/
│   ├── original_panel/                       Raw district- and region-level statistics (State Statistical Committee of Azerbaijan)
│   └── ...                                   Merged/processed analytical datasets
├── output/                                   Generated tables (.csv) -- 44 files
├── figures/                                  Generated figures (.png) -- Figures 2-4, S1-S2
├── graphical_abstract/
│   └── Graphical_Abstract.png                Manuscript graphical abstract
├── docs/
│   ├── CODEBOOK.md                           Analytical workflow, every model equation, script-to-output correspondence
│   ├── DATA_DESCRIPTION.md                   Data sources, provenance, and variable definitions
│   ├── REPRODUCIBILITY_CHECKLIST.md          Reproducibility checklist and internal-consistency verification against the manuscript
│   └── Replication_Guide.docx                Complete, step-by-step replication guide
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── .zenodo.json
├── LICENSE
├── requirements.txt
├── environment.yml
└── .gitignore
```

## Documentation

- **docs/CODEBOOK.md** -- analytical workflow, every model equation
  (DEA, SFA, Malmquist, second-stage, event-study), the control-frontier
  bootstrap adaptation, and the script-to-output correspondence
  (important -- read before extending this package)
- **docs/DATA_DESCRIPTION.md** -- data sources, provenance, and
  variable-level definitions for every model input
- **docs/REPRODUCIBILITY_CHECKLIST.md** -- reproducibility checklist and
  full internal-consistency verification against the manuscript,
  including a documented account of every analytical correction made
  during preparation of this package
- **docs/Replication_Guide.docx** -- complete replication guide

## Installation

```bash
conda env create -f environment.yml
conda activate agri-systems-repro
```

or

```bash
pip install -r requirements.txt
```

## Run

```bash
cd code
python run_all.py
```

This reproduces the complete analytical workflow: data preparation,
bootstrap DEA, SFA, second-stage analysis, bootstrap Malmquist
productivity decomposition, the interference-robust (control-frontier)
event-study estimation with clustered standard errors and
minimum-detectable-effect calculation, augmented-technology sensitivity
checks, robustness diagnostics, the additional review-driven robustness
checks, and all main-text and Supporting Information figures.

Expected total runtime is approximately 60-75 minutes, dominated by four
bootstrap procedures (see the docstring in `code/run_all.py` for a
per-script breakdown). Each bootstrap script prints its own progress to
the console and writes results incrementally, so a run can be monitored
or resumed script-by-script if interrupted.

`code/08a_make_figure1_framework.py` is retained in the repository but
is **not** part of the active pipeline and is not called by
`run_all.py`. It produces an analytical-framework diagram that is not
referenced anywhere in the current manuscript body; see
`docs/CODEBOOK.md` for the full note.

## Script-to-Output Correspondence

| Script | Produces |
|---|---|
| `00_data_prep.py` | Analytical panels; Table 1 sample-accounting assertion check (console) |
| `01_dea_bootstrap.py` | Table 2 (DEA summary), Figure 2 -- **core efficiency analysis** |
| `02_sfa_models.py` | Table 3 (translog SFA), SI Tables S1 (Cobb-Douglas), S2 (panel trend) |
| `03_second_stage.py` | Table 5 (determinants), SI Tables S3.1 (lagged), S5/S5b (static comparison), S6 (reduced form), S9 (outlier sensitivity) |
| `04_malmquist_bootstrap.py` | Table 4 (Malmquist decomposition), Figure 3, SI Figures S1-S2 -- **core productivity-growth analysis** |
| `05_event_study.py` | Legacy diagnostic only; not the source of any reported result |
| `06_labor_capital_land_robustness.py` | Table 9 (augmented-technology checks), SI Tables S10-S12 |
| `07_robustness_diagnostics.py` | SI Tables S4 (balanced panel), S7 (VIF), S8 (district ranking) |
| `08_make_figures.py` | Figures 2-4 (main text), Figures S1-S2 (SI) |
| `09_control_frontier_event_study.py` | Tables 6-8, Figure 4 -- **core causal-evaluation analysis** |
| `10_gpt_review_robustness.py` | SI Tables S15 (rank correlation), S16 (year-FE second stage), S17 (control-frontier diagnostics), S18 (Aghdara exclusion + wild cluster bootstrap) |
| `run_all.py` | All of the above, in dependency order |

## A Note on Script Execution Order and Dependencies

Scripts must be run in numerical order the first time, because several
depend on the output of earlier ones (`03`, `06`, `07`, `09`, and `10`
all read `output/table3_dea_bootstrap_full.csv`, produced by `01`; `10`
additionally reads outputs of `02`, `03`, and `09`). `run_all.py`
enforces this order automatically. Running an individual numbered
script directly is supported for regenerating a single part of the
pipeline, provided its own upstream dependencies have already been run
at least once.

## A Note on the Control-Frontier Design (Tables 6-8)

The manuscript's primary causal-evaluation specification
(`09_control_frontier_event_study.py`) scores every district's
technical efficiency against a reference technology built exclusively
from the 58 never-treated districts, so that no treated district's
post-2020 production choices can mechanically shift the frontier
against which any district -- treated or control -- is measured. This
is a direct methodological response to a concern (raised in an earlier
round of external review) that a difference-in-differences design
applied to *pooled*-frontier DEA scores can conflate a genuine treatment
effect with a purely mechanical artefact of frontier construction. See
`docs/CODEBOOK.md`, "Why a Separate Control-Frontier Script?", for the
full account, including how the bootstrap bias-correction procedure was
adapted for this asymmetric (58-district reference vs. 67-district
scored) construction.

## Citation

Please cite both the published article (once available) and this
archived repository. Citation metadata are provided in `CITATION.cff`
and `.zenodo.json`.

## License

MIT License (code and derived/processed data in this repository). See
`docs/DATA_DESCRIPTION.md` for the provenance of the underlying
district-level statistics and their own terms of use.

## Contact

**Halil Tosun**

Department of Animal Science, School of Agricultural and Food Sciences,
ADA University, Baku, Azerbaijan

ORCID: https://orcid.org/0000-0001-5117-0390

Email: halilibrahimtosun@gmail.com

**Matteo Vittuari**

Department of Agricultural and Food Sciences (DISTAL), University of
Bologna, Italy

School of Agricultural and Food Sciences, ADA
University, Baku, Azerbaijan

ORCID: https://orcid.org/0000-0003-4327-1575

Email: mvittuari@ada.edu.az

**Zenodo DOI:** https://doi.org/10.5281/zenodo.21380280

**Version:** 2.0.0
