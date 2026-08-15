# CODEBOOK

## Analytical Workflow

This package estimates technical efficiency and productivity growth in
Azerbaijan's dairy sector (2000-2024), then evaluates whether the
post-2020 territorial reintegration of nine districts is associated with
a measurable, robustly identified change in dairy-sector technical
efficiency, within a single internally consistent panel.

| Script | Description | Produces |
|---|---|---|
| `00_data_prep.py` | Loads and merges raw district- and region-level statistics; builds the DEA analysis sample; asserts the Table 1 sample-accounting subset relationship (see "Table 1 Sample-Accounting Check" below) | Analytical panels; console assertion check |
| `01_dea_bootstrap.py` | **Core efficiency analysis.** Input-oriented DEA (VRS and CRS) for each annual cross-section, with Simar-Wilson (1998, 2000) smoothed bootstrap bias correction, B=200; includes the boundary-case degenerate-draw filter (see below) | Table 2, Figure 2 |
| `02_sfa_models.py` | Translog and Cobb-Douglas stochastic frontier models (normal-half-normal), and a panel specification with a linear time trend | Table 3, SI S1, S2 |
| `03_second_stage.py` | Simar-Wilson (2007) truncated-normal second-stage regression of bootstrap-corrected DEA efficiency on district-level covariates; lagged, reduced-form, and outlier-treatment sensitivity variants | Table 5, SI S3.1, S5, S6, S9 |
| `04_malmquist_bootstrap.py` | **Core productivity-growth analysis.** Bootstrap Malmquist index (Simar and Wilson, 1999), B=300, decomposing cumulative TFP change into technological change (TC) and technical efficiency change (TEC) for 13 economic regions | Table 4, Figure 3, SI Figures S1-S2 |
| `05_event_study.py` | Legacy pooled-frontier event-study implementation, superseded by `09`; retained in the repository for internal comparison only and not the source of any value reported in the manuscript | -- |
| `06_labor_capital_land_robustness.py` | Labor-, capital-, and land-augmented DEA/SFA sensitivity checks on independently sourced, restricted subsamples | Table 9, SI S10-S12 |
| `07_robustness_diagnostics.py` | Balanced-panel comparison, variance inflation factors, district-level efficiency ranking | SI S4, S7, S8 |
| `08_make_figures.py` | Generates all main-text and Supporting Information figures at the exact plotting parameters (300 DPI) used in the submitted manuscript; the single authoritative figure-generation script for this package | Figures 2-4, S1-S2 |
| `08a_make_figure1_framework.py` | Generates a conceptual analytical-framework diagram. **Not called by `run_all.py`** -- this figure is not referenced anywhere in the current manuscript body; retained for reference only | (unused) |
| `09_control_frontier_event_study.py` | **Core causal-evaluation analysis.** Control-only-frontier bootstrap DEA construction (see "Why a Separate Control-Frontier Script?" below); pooled two-way fixed-effects DiD and dynamic event-study estimates under classical, HC1-robust, and district-clustered standard errors; fractional-logit (quasi-binomial) specification; minimum-detectable-effect and power calculations | Tables 6-8, Figure 4 |
| `10_gpt_review_robustness.py` | Four additional robustness/diagnostic analyses added following external methodological review: DEA-SFA rank correlation; second-stage regression with year fixed effects; control-frontier bootstrap and raw-theta-truncation diagnostics; Aghdara-exclusion sensitivity check and wild cluster (Rademacher) bootstrap on the primary causal estimate | SI S15-S18 |
| `run_all.py` | Runs scripts `00`-`10` in dependency order, including `05` as a legacy diagnostic comparison (`08a` is not called; see below) | All tables and figures |

## Why a Separate Control-Frontier Script? (Tables 6-8)

Script `01_dea_bootstrap.py` constructs, for each year, a reference
technology from **all 67 districts pooled together**, which is
appropriate for the descriptive efficiency levels reported in Table 2
but creates a specific problem for causal evaluation: because the nine
treated (reintegrated) districts' post-2020 production choices are part
of the pooled reference technology, their inclusion can, in principle,
mechanically shift the measured efficiency of every district -- treated
and untreated alike -- scored against that frontier in a given year.
This is a violation of the Stable Unit Treatment Value Assumption
(SUTVA) that a difference-in-differences design applied directly to
pooled-frontier scores cannot rule out.

Script `09_control_frontier_event_study.py` instead constructs, for
each year, the reference technology **exclusively from the 58
never-treated districts**, then scores every district (treated and
control) against this externally defined technology. By construction,
no treated district's post-2020 production can move this frontier. The
standard Simar-Wilson bootstrap (which perturbs and re-scores the same
set of units) does not apply directly to this asymmetric
(reference-set-differs-from-scored-set) construction and is adapted as
follows: in each year, the bootstrap perturbs only the 58 control
districts' own input-output data, and every district -- treated and
control alike -- is re-scored against each perturbed control-only
reference set. This is the paper's **primary** causal-evaluation
specification; the pooled-frontier specification of `01` is reported
alongside it only as a secondary comparison (Table 6).

## Table 1 Sample-Accounting Check

`00_data_prep.py` includes a runtime `assert` verifying that the
1,523-observation complete-case production panel used throughout the
paper is *exactly* the set of district-years for which fodder-crop sown
area is recorded (the more restrictive of the three production
variables), and that every one of those 1,523 district-years also has
non-missing milk production and cattle stock data, with zero exceptions.
This directly explains why Table 1 reports N=1,648 for milk production
and cattle stock individually but N=1,523 for the analytical sample used
in every model in the paper (see Results, Section 3.1, and the table
notes). If this assertion ever fails on a data update, `00_data_prep.py`
will raise an error rather than silently propagate an unexplained
sample-accounting discrepancy downstream.

## Core Equations

| Quantity | Equation | Source |
|---|---|---|
| Input-oriented DEA (VRS/CRS) | θ_i = min θ s.t. θx_i - Xλ ≥ 0, Yλ ≥ y_i, (1'λ=1 for VRS), λ≥0 | Banker, Charnes and Cooper, 1984 (VRS); Charnes, Cooper and Rhodes, 1978 (CRS) |
| Bootstrap bias correction | θ̃_i = 2θ̂_i - (1/B)Σ_b θ*_{i,b}, truncated to [0,1] | Simar and Wilson, 1998, 2000 |
| Control-frontier reference technology | θ̂_{i,control} = min θ s.t. θx_i - X_c λ ≥ 0, Y_c λ ≥ y_i, 1'λ=1, λ≥0 (X_c, Y_c = 58 control districts only) | This study, adapting Simar and Wilson (1998, 2000) |
| Translog SFA | ln y_i = β_0 + Σ_k β_k ln x_ki + (1/2)Σ_k Σ_l β_kl ln x_ki ln x_li + v_i - u_i | Aigner, Lovell and Schmidt, 1977 |
| Technical efficiency prediction | TE_i = E[exp(-u_i) \| ε_i] | Jondrow et al., 1982; Battese and Coelli, 1988 |
| Bootstrap Malmquist index | M_i = [ (D^t(x_{t+1},y_{t+1})/D^t(x_t,y_t)) × (D^{t+1}(x_{t+1},y_{t+1})/D^{t+1}(x_t,y_t)) ]^(1/2) = TEC × TC | Färe et al., 1994; bootstrap: Simar and Wilson, 1999 |
| Second-stage truncated regression | θ̃_i = z_i'β + ε_i, ε_i ~ N(0,σ_ε²) truncated at (1-z_i'β) | Simar and Wilson, 2007 |
| Pooled DiD | θ̃_{it} = α_i + λ_t + δ(Treated_i × Post_t) + ε_it | Standard two-way fixed effects |
| Event study | θ̃_{it} = α_i + λ_t + Σ_{k=-10,k≠-1}^{4} β_k(Treated_i × 1[t=2020+k]) + ε_it | Standard dynamic DiD |
| Minimum detectable effect | MDE = (z_0.975 + z_0.80) × SE ≈ 2.80 × SE | Bloom, 1995; List, Sadoff and Wagner, 2011 |
| Wild cluster bootstrap | Rademacher weights on the null-restricted model, B=999 | Cameron, Gelbach and Miller, 2008; Cameron and Miller, 2015 |

Full bibliographic details, including DOIs where available, are
provided in the accompanying manuscript's References section.

## Boundary-Case Bootstrap Replicates (Degenerate-Draw Filter)

For district-years whose raw score lies exactly on the frontier
(θ̂_i = 1), the smoothed-reflection bootstrap can occasionally perturb
the pseudo-reference technology into a configuration for which the
input-oriented linear program has no economically meaningful solution
near 1; the solver then returns a numerically "successful" but
degenerate value (observed: exactly 1e6, an unbounded-variable ceiling
artefact) rather than failing outright. `01_dea_bootstrap.py` and
`09_control_frontier_event_study.py` both filter any bootstrap replicate
exceeding `DEGENERATE_THETA_CEILING = 5` and exclude it from that unit's
bias-correction mean via `nanmean`; if every one of a unit's B replicates
is degenerate (rare, and more common under the smaller 58-district
control-only reference set than under the full 67-district pooled
reference set), the raw score is retained without correction rather than
reported as missing. Full per-observation diagnostics are in
`output/table3_dea_bootstrap_full.csv` (columns
`vrs_bootstrap_degenerate_draws`) and SI Table S17
(`output/s17a_control_frontier_diagnostics_by_year.csv`,
`output/s17b_raw_theta_gt1_by_group.csv`).

## Determinism

`01_dea_bootstrap.py`, `04_malmquist_bootstrap.py`, and
`09_control_frontier_event_study.py` each use a fixed random seed
(`np.random.default_rng(seed)` with an explicit integer seed passed at
each call site) for every bootstrap procedure. Given the same input data
and the same NumPy/SciPy versions, every script in this package produces
output that reproduces to the precision reported in the manuscript on
every run and on every machine; see
`docs/REPRODUCIBILITY_CHECKLIST.md` for an exact, independently
re-verified account, including a pixel-identical (MD5 checksum) match
for all five figures across two independent full-pipeline runs.

## Optimization / Estimation Methods

- **DEA linear programs**: `scipy.optimize.linprog` (HiGHS method).
- **Maximum-likelihood estimation** (SFA, truncated regression,
  fractional logit): Nelder-Mead followed by BFGS refinement, implemented
  directly against SciPy (no higher-level statistical package dependency,
  so every estimator's likelihood, gradient handling, and standard-error
  calculation is fully transparent in this source code).
- **Cluster-robust and wild-cluster-bootstrap standard errors**: direct
  sandwich-estimator implementation (`09_control_frontier_event_study.py`,
  `10_gpt_review_robustness.py`); not dependent on an external
  econometrics package.
