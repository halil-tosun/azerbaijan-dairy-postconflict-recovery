# Understanding Post-conflict Agricultural Systems Recovery

## Open Science Replication Package

This repository contains the complete replication package accompanying the manuscript **"Understanding Post-conflict Agricultural Systems Recovery: An Integrated Efficiency, Productivity and Policy Analysis of Azerbaijan's Dairy Sector"** prepared for submission to **Agricultural Systems (Elsevier)**.

---

## Repository Overview

This repository follows open science and computational reproducibility principles and includes:

- Complete Python source code
- Raw and analytical datasets
- Reproducible computational workflow
- Comprehensive documentation
- Software environment specifications
- Replication guide

---

## Repository Structure

```text
Azerbaijan_Dairy_Replication/
├── code/
├── data/
├── docs/
├── manuscript/
├── output/
├── README.md
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── environment.yml
└── .gitignore
```

## Documentation

- **docs/CODEBOOK.md** – Analytical workflow and variable descriptions
- **docs/DATA_DESCRIPTION.md** – Data sources and dataset structure
- **docs/REPRODUCIBILITY_CHECKLIST.md** – Reproducibility checklist
- **docs/Replication_Guide.docx** – Complete replication guide

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
python code/run_all.py
```

This reproduces the complete analytical workflow including data preparation, DEA, SFA, second-stage analysis, Malmquist productivity decomposition, event-study estimation, robustness analyses, and figure generation.

## Citation

Please cite both the published article and the archived GitHub/Zenodo repository. Citation metadata are provided in `CITATION.cff`.

## License

MIT License.

## Contact

**Halil Tosun, Ph.D.**

ORCID: https://orcid.org/0000-0001-5117-0390

Email: halilibrahimtosun@gmail.com

**Version:** 1.0.0

**Zenodo DOI:** To be assigned after public release.
