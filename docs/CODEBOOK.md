# CODEBOOK

## Replication Package

**Understanding Post-conflict Agricultural Systems Recovery: An Integrated Efficiency, Productivity and Policy Analysis of Azerbaijan's Dairy Sector**

**Author:** Halil Tosun, Ph.D.

---

## Purpose

This document describes the role of each script included in the replication package, its required inputs, and its expected outputs. It serves as a roadmap for reproducing all empirical analyses reported in the manuscript.

---

## Workflow

```text
Raw Data
    |
    v
00_data_prep.py
    |
    v
Processed Dataset
    |
    +--> 01_dea_bootstrap.py
    |         |
    |         v
    |    03_second_stage.py
    |
    +--> 02_sfa_models.py
    |
    +--> 04_malmquist_bootstrap.py
    |
    +--> 05_event_study.py
    |
    +--> 06_labor_capital_land_robustness.py
    |
    +--> 07_robustness_diagnostics.py
    |
    +--> 08a_make_figure1_framework.py
    |
    +--> 08_make_figures.py
    |
    v
run_all.py
```

## Script Descriptions

### 00_data_prep.py
Prepares the raw agricultural datasets and creates the processed analytical panel.

### 01_dea_bootstrap.py
Estimates bootstrap Data Envelopment Analysis (DEA) efficiency scores.

### 02_sfa_models.py
Estimates Stochastic Frontier Analysis (SFA) models.

### 03_second_stage.py
Implements Simar–Wilson second-stage truncated regression.

### 04_malmquist_bootstrap.py
Calculates bootstrap Malmquist productivity indices and decomposition.

### 05_event_study.py
Estimates dynamic treatment effects using an event-study framework.

### 06_labor_capital_land_robustness.py
Performs robustness analyses using alternative production factor specifications.

### 07_robustness_diagnostics.py
Runs additional robustness and diagnostic analyses.

### 08a_make_figure1_framework.py
Generates the conceptual framework figure.

### 08_make_figures.py
Creates the remaining manuscript figures.

### run_all.py
Executes the complete replication workflow.

Run:

```bash
python run_all.py
```

to reproduce the complete set of analyses, tables, and figures.

---

## Expected Outputs

The workflow generates:

- Processed datasets
- DEA and SFA efficiency estimates
- Malmquist productivity indices
- Regression outputs
- Manuscript tables
- Manuscript figures

Outputs are saved under the `output/` directory.

---

## Reproducibility

The repository is organized to maximize transparency and computational reproducibility. All analyses can be reproduced from the supplied data and source code.
