#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
TREE_RESULTS=${ROOT}/results/01_results_phylo_tree
RESULTS=${ROOT}/results/02_results_div_times
LOGS=${ROOT}/logs/02_logs_div_times

TREE_FILE=${TREE_RESULTS}/phylogeny_rooted.treefile



### create results/log folder(s)
mkdir -p "${RESULTS}"
mkdir -p "${LOGS}"



### setup logging
JOBTAG="${SLURM_JOB_ID}_div_times"
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
echo "Date: $(date)"
echo "Tree file: ${TREE_FILE}"
echo "Results directory: ${RESULTS}"



### sanity checks
if [ ! -f "${TREE_FILE}" ]; then
    echo "ERROR: Tree file not found: ${TREE_FILE}"
    echo "Run 01_phylogenetic_tree.sh first"
    exit 1
fi



### ===========================================================================
### STEP 1: Prepare taxonomy file for PhyloRank
### ===========================================================================
echo ""
echo "### STEP 1: Preparing taxonomy file ###"
echo ""

conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

python3 << EOF
from pathlib import Path
import re

tree_file = Path("${TREE_FILE}")
tax_file = Path("${RESULTS}/taxonomy.tsv")

with open(tree_file) as f:
    tree_str = f.read()

tips = re.findall(r'([A-Za-z0-9_.-]+):', tree_str)
print(f"Found {len(tips)} tips in tree")

with open(tax_file, "w") as out:
    for tip in tips:
        if tip.startswith("Pseflo"):
            taxonomy = "d__Fungi;p__Ascomycota;c__Dothideomycetes;o__Gloniales;f__Gloniaceae;g__Pseudocenococcum;s__Pseudocenococcum floridanum"
        else:
            taxonomy = "d__Fungi;p__Ascomycota;c__Dothideomycetes;o__Gloniales;f__Gloniaceae;g__Cenococcum;s__Cenococcum geophilum"
        out.write(f"{tip}\t{taxonomy}\n")

print(f"Taxonomy file written to: {tax_file}")
EOF

echo ""
echo "STEP 1 complete."

conda deactivate



### ===========================================================================
### STEP 2: Decorate tree with PhyloRank (taxonomy + consistency, no RED)
### ===========================================================================
echo ""
echo "### STEP 2: Decorating tree with PhyloRank ###"
echo ""

conda activate HTT_phylorank

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

phylorank decorate \
    "${TREE_FILE}" \
    "${RESULTS}/taxonomy.tsv" \
    "${RESULTS}/decorated_tree" \
    --skip_rd_refine

if [ $? -ne 0 ]; then
    echo "ERROR: phylorank decorate failed"
    exit 1
fi

echo ""
echo "Decorate complete."
echo "Output files:"
ls -la "${RESULTS}"/decorated_tree*

conda deactivate



### ===========================================================================
### STEP 3: Manually calculate RED values from the rooted tree
### ===========================================================================
echo ""
echo "### STEP 3: Calculating RED values manually (Parks et al. formula) ###"
echo ""

conda activate HTT_python_ete3

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

python3 << EOF
from ete3 import Tree
from pathlib import Path
import math

tree_file = Path("${TREE_FILE}")
out_file  = Path("${RESULTS}/red_values.tsv")

print(f"Reading tree from: {tree_file}")
t = Tree(str(tree_file), format=1)

# Ensure branch lengths are defined
for n in t.traverse():
    if n.dist is None:
        n.dist = 0.0

red = {}

def mean_parent_to_tips(parent, child):
    """Mean distance from parent to all descendant tips of child."""
    tips = child.get_leaves()
    if not tips:
        return 0.0
    dists = [parent.get_distance(leaf) for leaf in tips]
    return sum(dists) / len(dists)

root = t.get_tree_root()
red[root] = 0.0  # RED(root) = 0

# Preorder traversal: ensure parent RED computed before children
for parent in root.traverse("preorder"):
    p_red = red[parent]
    for child in parent.children:
        d = child.dist
        u = mean_parent_to_tips(parent, child)
        if u <= 0 or (not math.isfinite(u)):
            # Fallback: propagate parent RED if degenerate
            red_child = p_red
        else:
            red_child = p_red + (d / u) * (1.0 - p_red)
        red[child] = red_child

# Name tips with their actual names and give unique names to all internal nodes
rows = []
idx = 1
for n in t.traverse():
    if n.is_leaf():
        name = n.name
    else:
        name = f"internal_{idx}"
        idx += 1
    rows.append((name, red[n], n.is_leaf()))

out_file.parent.mkdir(parents=True, exist_ok=True)
with out_file.open("w") as out:
    out.write("node\tred\tis_tip\n")
    for name, val, is_tip in rows:
        out.write(f"{name}\t{val:.6f}\t{int(is_tip)}\n")

print(f"Wrote RED values to: {out_file}")

# Quick sanity check
tip_reds = [val for n, val, is_tip in rows if is_tip]
if tip_reds:
    print(f"Tip RED range: {min(tip_reds):.4f} - {max(tip_reds):.4f}")
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: manual RED calculation failed"
    exit 1
fi

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
echo "  Taxonomy file:            ${RESULTS}/taxonomy.tsv"
echo "  Decorated tree:           ${RESULTS}/decorated_tree"
echo "  Statistics table:         ${RESULTS}/decorated_tree-table"
echo "  Summary table:            ${RESULTS}/decorated_tree-summary"
echo "  Taxonomy consistency:     ${RESULTS}/decorated_tree-taxonomy"
echo "  RED values (manual):      ${RESULTS}/red_values.tsv"