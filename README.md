# Horizontal transfer of transposable elements across lineages of an ectomycorrhizal fungal species complex


This repository contains all scripts, software environment definitions, and directory structure required to reproduce the analyses presented in:

> **...**  

---

## Abstract

Within closely related, putatively non-recombining species, the role of horizontal transfer of transposable elements (HTT) remains poorly understood, while it is increasingly recognized as a contributor to genome evolution in eukaryotes. Here, we investigate HTT in a clonal fungal species complex characterized by pronounced genome size variation of unknown origin. Using a phylogeny-aware framework applied to 22 genome assemblies, we identified candidate HTT events by comparing synonymous substitution rates of transposable element protein domains to lineage-specific divergence thresholds estimated from conserved single-copy orthologs. To conservatively account for shared ancestry and repeated transfers, candidate HTTs were clustered into independent transfer events and reduced to a minimal set of parsimonious HTT events across the species tree.

We detected thousands of HTT candidates spanning multiple TE superfamilies, with strong enrichment in specific clades and element types. Synteny-aware analysis of genomic flanking regions confirmed that HTT-associated TEs lack conserved sequence context beyond element boundaries, consistent with independent insertion rather than vertical inheritance or introgression. HTT frequency varied markedly among lineages, suggesting that transfer susceptibility is shaped by clade-specific factors rather than genome or mobilome size.

Together, our results demonstrate that horizontal transfer of transposable elements occurs even among closely related, clonal fungal lineages and represents a source of genetic exchange and mobilome turnover in the absence of sexual recombination.

---

## Repository Structure

```
2026_Cenge_HTT/
├── HPC/                         # High-performance computing (cluster)
│   ├── envs/                    # Conda environment definitions
│   └── scripts/                 # Shell / batch scripts for HPC execution
│
├── R/                           # Local R-based analyses
│   ├── data/                    # Inputs for R scripts (derived from HPC results)
│   ├── results/                 # Tables and figures generated in R
│   ├── scripts/                 # R analysis scripts
│   └── HTT.Rproj                # RStudio project file
│
└── README.md
```

---

## Workflow Overview

The analysis consists of two components:

1. **HPC-based analyses**  
   Computationally intensive steps (phylogeny inference, divergence estimation, TE annotation, HTT detection, clustering, and flanking region analysis) are executed on a computing cluster.

2. **Local R-based analyses**  
   Downstream statistical analyses and figure generation are performed locally using R.

---

## HPC: Cluster-Based Analyses

### Setup

From the repository root:

```bash
cd HPC
bash scripts/00_setup.sh
```

This will:
- Create all required Conda environments
- Initialize directory structure for logs and results

### Execution

HPC scripts are located in:

```
HPC/scripts/
```

Scripts are intended for submission to a job scheduler (e.g. SLURM). Scheduler directives may need to be adapted to your system.

---

## R: Local Analyses

Open the RStudio project:

```
R/HTT.Rproj
```

Run scripts in `R/scripts/` in numerical order.  
Inputs are read from `R/data/` and outputs written to `R/results/`.

---

## Required External Databases

Several analyses in this pipeline rely on external databases that are not distributed with this repository due to size, update frequency, and licensing considerations. Users must ensure that the appropriate databases are installed and accessible prior to running the full pipeline.

### Pfam

The Pfam protein domain database is used for domain-based annotation of transposable element proteins.

Users should download an appropriate Pfam release (e.g. Pfam-A HMMs) from:
https://ftp.ebi.ac.uk/pub/databases/Pfam/

After download, the database must be prepared for use with HMMER (e.g. using `hmmpress`).

### CDD (Conserved Domain Database)

The NCBI Conserved Domain Database (CDD) is used for complementary protein domain annotation and classification.

CDD data can be obtained from:
https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml

Depending on the analysis step, users may need either the CDD PSSMs or a locally formatted database compatible with RPS-BLAST or related tools.

### EarlGrey Databases

Transposable element annotation using **EarlGrey** requires additional reference databases to be installed and properly configured. These typically include repeat libraries and auxiliary databases required by the underlying tools invoked by EarlGrey.

Users should follow the official EarlGrey documentation for database installation and setup:
https://github.com/TobyBaril/EarlGrey

All required databases must be downloaded and indexed prior to running EarlGrey-based annotation steps. Paths to these databases are specified within the corresponding HPC scripts and may need to be adapted to local system layouts.

### Database Location and Configuration

The pipeline assumes that Pfam, CDD, and EarlGrey-related databases are available locally and accessible to the relevant HPC scripts. Database paths are defined within the scripts or configuration sections and should be adjusted as needed.

A typical local directory structure may resemble:

```
databases/
├── pfam/
├── cdd/
└── earlgrey/
```

databases/
├── pfam/
├── cdd/
└── earlgrey/

Users running the pipeline on shared HPC systems may alternatively point to centrally maintained database installations.

### Reproducibility Note

Exact database versions and releases used in the original analysis are documented in the associated Zenodo archive and/or supplementary materials. For full reproducibility, users are encouraged to use the same database versions where possible.

---

## Data and Results Availability

Due to size constraints, large input data and intermediate/results files are not stored directly in this GitHub repository. These data are archived on Zenodo in two complementary records corresponding to the two components of the workflow.

### HPC Data and Results

All data required to run or inspect the cluster-based analyses (including intermediate files and results generated on the HPC system) are archived on Zenodo:

🔗 https://doi.org/XXXX.XXXX/zenodo.HPC_DOI

This archive corresponds to the directory structure expected by the scripts in `HPC/`. After downloading, the contents should be placed into the `HPC/` directory as indicated in the archive README.

### R Data and Results

All input data and results used for downstream statistical analyses and figure generation in R are archived separately on Zenodo:

🔗 https://doi.org/XXXX.XXXX/zenodo.R_DOI

This archive corresponds to the `R/data/` and `R/results/` directories. After downloading, users should place the archived contents into the respective directories before running the R scripts.

### Relationship Between GitHub and Zenodo Archives

The GitHub repository contains all scripts, software environment definitions, and directory structure required to reproduce the analyses. The Zenodo archives provide the large data files required to either:

1. Reproduce all figures and tables directly from archived results, or  
2. Inspect, validate, or extend the analyses without re-running the full HPC pipeline.

Users who wish to fully recompute the analyses from raw inputs may instead follow the workflow described in the repository and obtain genome assemblies as described below.

---

## Reproducibility

- All software dependencies are defined via Conda environments
- Analysis steps are explicitly ordered and numbered
- Directory structure is created automatically
- Scripts are safe to rerun

---

## Citation

If you use this repository, please cite:

*...*
