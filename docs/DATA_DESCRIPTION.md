# DATA_DESCRIPTION

## Overview

This document describes the datasets included in the replication package for:

**Technical Efficiency, Productivity Growth, and the Effects of Territorial Reintegration: Evidence from Azerbaijan's Dairy Sector, 2000–2024**

Author: **Halil Tosun, Ph.D.**

---

## Data Sources

The empirical analyses are based on official agricultural statistics for Azerbaijan compiled into panel datasets for efficiency, productivity, and policy evaluation.

The repository separates raw source files from analytical panel datasets.

---

## Repository Structure

```text
data/
├── official_stats_raw/
│   ├── CSV and XLS files obtained from official statistical sources
│   └── Used as the original source data
│
└── original_panel/
    ├── azerbaijan_dairy_panel_long.csv
    ├── azerbaijan_dairy_panel_NATIONAL.csv
    ├── azerbaijan_dairy_panel_WIDE_districts.csv
    ├── azerbaijan_dairy_panel_WIDE_economic_regions.csv
    ├── azerbaijan_DISTRICTS_final_efficiency_outputs.csv
    └── azerbaijan_national_milk_yield_per_cow-1.csv
```

---

## Raw Data

The **official_stats_raw/** directory contains the original statistical tables in both CSV and XLS formats. These files constitute the primary data sources used to construct the analytical datasets.

No modifications should be made to these files.

---

## Analytical Data

The **original_panel/** directory contains harmonized datasets used throughout the empirical analyses.

These datasets provide the basis for:

- Data preparation
- DEA estimation
- SFA estimation
- Second-stage analysis
- Malmquist productivity decomposition
- Event-study estimation
- Robustness analyses
- Figure generation

---

## Temporal Coverage

2000–2024

---

## Geographic Coverage

Republic of Azerbaijan, including district, economic-region, and national aggregation levels where applicable.

---

## Reproducibility Notes

- Raw data are preserved separately from analytical datasets.
- All empirical results can be reproduced using the supplied Python scripts.
- Analytical outputs are generated automatically by running the replication workflow.

---

## Citation

If these data are used, please cite both the associated journal article and the archived GitHub/Zenodo repository.
