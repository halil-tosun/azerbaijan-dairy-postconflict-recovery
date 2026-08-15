# CHANGELOG

All notable changes to this replication package will be documented in this file.

The format is inspired by *Keep a Changelog* and follows semantic versioning where appropriate.

---

## Version 2.0.0


### Added
- `code/09_control_frontier_event_study.py`: interference-robust ("control-only
  frontier") event-study design addressing a Stable Unit Treatment Value
  Assumption (SUTVA) concern in the original pooled-frontier design; adds
  district-clustered standard errors, a fractional-logit (quasi-binomial)
  specification, and minimum-detectable-effect / power calculations.
  Produces Tables 6-8 and Figure 4 of the revised manuscript.
- `code/10_gpt_review_robustness.py`: four additional robustness/diagnostic
  analyses added in response to external methodological review — DEA-SFA
  rank correlation (S15), second-stage regression with year fixed effects
  to rule out a nominal-price confound (S16), control-frontier bootstrap
  and raw-theta-truncation diagnostics (S17), and an Aghdara-exclusion
  sensitivity check plus a wild cluster (Rademacher) bootstrap on the
  primary causal estimate (S18).
- `code/08_make_figures.py`: consolidated, single authoritative figure-
  generation script (300 DPI) for all main-text and Supporting Information
  figures, replacing an earlier version whose output filenames and DPI did
  not match the manuscript.
- `graphical_abstract/Graphical_Abstract.png`: manuscript graphical abstract.
- `.zenodo.json`: Zenodo archival metadata (was previously missing from
  this repository).
- `code/00_data_prep.py`: runtime assertion verifying the Table 1
  sample-accounting subset relationship (every district-year with
  fodder-crop-area data also has non-missing milk production and cattle
  stock data); raises an error rather than silently propagating an
  unexplained N=1,648-vs-N=1,523 discrepancy if this relationship is
  ever violated by a future data update.
- `docs/REPRODUCIBILITY_CHECKLIST.md`: full internal-consistency
  verification table (every reported statistic in Tables 1-9 and S1-S18
  independently re-derived from a cleared output directory and compared
  against the manuscript), rewritten to this standard for the present
  release.

### Fixed
- `code/01_dea_bootstrap.py`: fixed a numerical defect in the Simar-Wilson
  bootstrap bias correction affecting district-years whose raw VRS score is
  exactly on the frontier (theta_hat = 1). The smoothed-reflection bootstrap
  could, for a fraction of replicates, return a numerically degenerate
  solution (observed: an unbounded-variable ceiling value) rather than a
  genuine failure; left unfiltered, this either silently dropped affected
  district-years (NaN) or corrupted their bias-correction mean. A
  `DEGENERATE_THETA_CEILING` filter now excludes degenerate replicates from
  the bias-correction mean, with a raw-score fallback for the rare case
  where all 200 replicates are degenerate. This resolved a previously
  unexplained "off-by-a-few" discrepancy between reported and reproduced N
  in several tables.
- `code/03_second_stage.py`: fixed a column-name collision in
  `build_second_stage_sample()` where the DEA-bootstrap output's own
  `cost_milk`/`profitability` columns were not dropped before merging in
  a (possibly differently outlier-treated) copy of the same columns,
  silently defeating the `outlier_treatment` argument and, on newer pandas
  versions, raising a `KeyError` that broke the one-command `run_all.py`
  workflow entirely.
- Table 4 (translog SFA): corrected a decimal-point transcription error in
  the manuscript (the cattle-stock coefficient is 1.4126, not 14.126; the
  underlying code output was always correct — only the manuscript text
  was wrong).
- Milk production cost unit corrected to AZN/centner throughout the
  manuscript (was inconsistently given as AZN/ton in one table).

### Changed
- `code/run_all.py`: updated to include scripts 09 and 10 in the pipeline
  sequence; `08a_make_figure1_framework.py` removed from the active pipeline
  (retained in the repository but unused — it produces a figure not
  referenced anywhere in the current manuscript).
- README.md, CITATION.cff, docs/CODEBOOK.md, docs/DATA_DESCRIPTION.md, and
  docs/REPRODUCIBILITY_CHECKLIST.md rewritten to reflect the current
  manuscript title and target journal (previously referenced a superseded
  title and a different target journal), and restructured to match a
  higher external documentation standard (script-to-output correspondence
  tables, an exact internal-consistency verification table for every
  reported statistic, and an explicit "known, documented analytical
  corrections" section).

### Reproducibility
- Full pipeline independently re-run from a cleared output directory this
  version; every main-text and Supporting Information table was verified
  to reproduce exactly (within reporting precision), and all five figures
  were verified pixel-identical (MD5 checksum) between independent runs.
- Expected runtime increased to ~60-75 minutes (previously ~30-40 minutes)
  due to the additional control-frontier and review-driven bootstrap
  procedures in scripts 09 and 10.

---

## Version 1.0.0 (Initial Public Release)

### Added
- Complete Python source code for all empirical analyses.
- Raw and analytical datasets.
- README.md with repository overview and usage instructions.
- CODEBOOK.md describing the analytical workflow.
- DATA_DESCRIPTION.md documenting data sources and structure.
- REPRODUCIBILITY_CHECKLIST.md.
- Replication_Guide.docx.
- CITATION.cff for software citation.
- LICENSE file.
- requirements.txt.
- environment.yml.
- .gitignore.

### Reproducibility
- One-command workflow via `run_all.py`.
- Computational environment documented.
- Repository prepared for GitHub release and Zenodo archiving.

### Notes
Zenodo DOI: https://doi.org/10.5281/zenodo.21380280

