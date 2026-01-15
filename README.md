# Horizontal transfer of transposable elements across lineages of an ectomycorrhizal fungal species complex

This repository contains scripts, software environment definitions, and directory structure required to reproduce the analyses presented in:

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

## Figure Reproducibility and Workflow Dependencies

All figures and tables presented in the associated publication can be reproduced by running the R scripts provided in this repository using the archived input data available on Zenodo.

Specifically:
- The R-based analyses (`R/scripts/`) generate all figures and summary tables used in the manuscript.
- These scripts can be run locally without access to an HPC system, provided that the required input data are placed in `R/data/` as described above.

The HPC-based pipeline is required **only** to regenerate the intermediate input files consumed by the R scripts (e.g. phylogenies, divergence estimates, transposable element annotations, and HTT candidate tables). Users who wish to inspect or modify the upstream analyses may rerun the HPC pipeline; otherwise, this step is not necessary for reproducing the published figures.

This separation allows full reproduction of all figures on standard desktop or laptop systems while retaining the ability to recompute all upstream analyses when needed.

---

## HPC: Cluster-Based Analyses

### Clone the repository on the HPC

```bash
git clone https://github.com/USERNAME/2026_Cenge_HTT.git
cd 2026_Cenge_HTT/HPC
```

### Setup

From `2026_Cenge_HTT/HPC` run:

```bash
bash scripts/00_setup.sh
```

This will:
- Create all required Conda environments from `HPC/envs/*.yml`
- Initialize directory structure under `HPC/` for `data/`, `databases/`, `logs/`, and `results/`


### Conda Environment Configuration (Important)

The Conda environment YAML files provided in `HPC/envs/` include a `prefix:` field that specifies the installation path of each environment.

The `prefix:` paths in these files are **placeholders** and must be adapted to the local system before creating the environments. Users may either:

1. **Edit the `prefix:` field** in each YAML file to point to a valid Conda environment directory on their system, or  
2. **Remove the `prefix:` field entirely**, allowing Conda (or mamba) to install the environment in its default location.

For example, the following entry in a YAML file:

```yaml
prefix: /path/to/conda/envs/HTT_hmmer
```

must be replaced with an appropriate local path or removed prior to environment creation.

Failure to update or remove the `prefix:` field will result in Conda attempting to install the environment into a non-existent directory.


### Execution

HPC scripts are located in:

```
HPC/scripts/
```

Scripts are intended for submission to a job scheduler (e.g. SLURM). Scheduler directives may need to be adapted to your system.

---

## R: Local Analyses

### Clone the repository locally

```bash
git clone https://github.com/nnamremmizxilef/2026_Cenge_HTT.git
cd 2026_Cenge_HTT/R
```

Open the RStudio project:

```
HTT.Rproj
```

Run scripts in `R/scripts/` in numerical order.  
Inputs are read from `R/data/` and outputs written to `R/results/`.

---

## Required External Databases

Several analyses in this pipeline rely on external databases that are not distributed with this repository due to size, update frequency, and licensing considerations. Users must ensure that the appropriate databases are installed and accessible prior to running the full pipeline.

### Pfam

Pfam protein domain database (Pfam-A HMMs):

https://ftp.ebi.ac.uk/pub/databases/Pfam/

After download, the database should be prepared for HMMER (e.g. `hmmpress`).

### CDD (Conserved Domain Database)

NCBI Conserved Domain Database:

https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml

Depending on the analysis step, users may need either the CDD PSSMs or a locally formatted database compatible with RPS-BLAST or related tools.

### EarlGrey Databases

Transposable element annotation using **EarlGrey** requires additional reference databases (repeat libraries and auxiliary databases). Follow the official EarlGrey documentation:

https://github.com/TobyBaril/EarlGrey

### BUSCO databases

BUSCO lineage datasets must be installed locally. The `busco_downloads/` directory should contain the lineage datasets referenced in the manuscript (directory names matching the lineages used in the analyses).

### Database Location and Configuration

The pipeline assumes that Pfam, CDD, EarlGrey-related databases, and BUSCO lineage datasets are available locally and accessible to the relevant HPC scripts. Database paths are defined within scripts or configuration sections and should be adjusted as needed.

A typical local directory structure may resemble:

```
databases/
├── pfam/
├── cdd/
├── earlgrey/
└── busco_downloads/
```

Users running the pipeline on shared HPC systems may alternatively point to centrally maintained database installations by updating the corresponding path variables in the scripts.

---

## Data and Results Availability (Zenodo)

Large input data and intermediate/results files are not stored directly in this GitHub repository. These data are archived on Zenodo in two records corresponding to the two components of the workflow.

### HPC Data (`data_hpc.zip`)

🔗 https://doi.org/XXXX.XXXX/zenodo.HPC_DOI

After downloading and extracting `data_hpc.zip`, place its contents into:

```
2026_Cenge_HTT/
└── HPC/
    └── data/
        └── (contents of data_hpc.zip)
```

### R Data (`data_R.zip`)

🔗 https://doi.org/XXXX.XXXX/zenodo.R_DOI

After downloading and extracting `data_R.zip`, place its contents into:

```
2026_Cenge_HTT/
└── R/
    └── data/
        └── (contents of data_R.zip)
```

If the R archive also provides precomputed outputs, place them into `R/results/` as instructed in the Zenodo record.

---

## Reproducibility

- All software dependencies are defined via Conda environments
- Analysis steps are explicitly ordered and numbered
- Directory structures are created automatically
- Scripts are safe to rerun

---

## Citation

If you use this repository, please cite:

*...*
