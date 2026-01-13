#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
RESULTS=${ROOT}/results/01_results_phylo_tree
LOGS=${ROOT}/logs/01_logs_phylo_tree

### define outgroup for tree rooting
OUTGROUP="PsefloM405_1_AssemblyScaffolds_2024-02-18"



### setup logging
JOBTAG="${SLURM_JOB_ID}_phylo_tree_final"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"



### environment setup
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export LANGUAGE=C

eval "$(conda shell.bash hook)"



### define start time
START_TIME=$(date +%s)



### info
echo "Running on HYPERION - Date: $(date)"
echo "Results directory: ${RESULTS}"
echo "Outgroup: ${OUTGROUP}"



### navigate to results directory
cd "${RESULTS}"

### sanity check
if [ ! -f "concatenated.faa" ]; then
    echo "ERROR: concatenated.faa not found in ${RESULTS}"
    echo "Run the full phylogenetic pipeline first"
    exit 1
fi

echo "Found concatenated alignment: $(wc -c < concatenated.faa) bytes"



### ===========================================================================
### STEP 7: Build phylogenetic tree with IQ-TREE
### ===========================================================================
echo ""
echo "### STEP 7: Building phylogenetic tree with IQ-TREE ###"
echo ""

conda activate HTT_iqtree

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "IQ-TREE version: $(iqtree2 --version 2>&1 | head -1)"

iqtree2 \
    -s concatenated.faa \
    -B 1000 \
    -m MFP \
    -mset LG,WAG,JTT \
    -mrate G4,R4 \
    -T AUTO \
    --prefix phylogeny

echo ""
echo "STEP 7 complete."

conda deactivate



### ===========================================================================
### STEP 8: Root tree with outgroup
### ===========================================================================
echo ""
echo "### STEP 8: Rooting tree with outgroup: ${OUTGROUP} ###"
echo ""

conda activate HTT_r

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "R version: $(R --version | head -1)"

Rscript --vanilla << EOF
library(ape)

tree <- read.tree("phylogeny.treefile")
outgroup <- "${OUTGROUP}"

# Check if outgroup exists in tree
if (!outgroup %in% tree$tip.label) {
    cat("ERROR: Outgroup not found in tree\n")
    cat("Looking for:", outgroup, "\n")
    cat("Available tips:\n")
    cat(paste(tree$tip.label, collapse = "\n"), "\n")
    quit(status = 1)
}

cat("Rooting with outgroup:", outgroup, "\n")
tree_rooted <- root(tree, outgroup = outgroup, resolve.root = TRUE)
write.tree(tree_rooted, "phylogeny_rooted.treefile")
cat("Rooted tree saved to: phylogeny_rooted.treefile\n")
EOF

echo ""
echo "STEP 8 complete."

conda deactivate



### ===========================================================================
### FINISH
### ===========================================================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "Total runtime: $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Tree file (unrooted): ${RESULTS}/phylogeny.treefile"
echo "  Tree file (rooted):   ${RESULTS}/phylogeny_rooted.treefile"
echo "  IQ-TREE log:          ${RESULTS}/phylogeny.iqtree"
echo ""