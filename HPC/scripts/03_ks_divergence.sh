#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
TREE_RESULTS=${ROOT}/results/01_results_phylo_tree
RESULTS=${ROOT}/results/03_results_ks_divergence
LOGS=${ROOT}/logs/03_logs_ks_divergence
BUSCO_RUNS_DIR=${TREE_RESULTS}/busco_runs

# input files
TREE_FILE=${TREE_RESULTS}/phylogeny_rooted.treefile
BUSCO_SEQUENCES_DIR=${TREE_RESULTS}/sequences

### define outgroup (same as phylogenetic analysis)
OUTGROUP="PsefloM405_1_AssemblyScaffolds_2024-02-18"



### create results/log folder(s)
mkdir -p "${RESULTS}"/{node_analyses,alignments,codon_alignments,ks_calculations}
mkdir -p "${LOGS}"



### setup logging
JOBTAG="${SLURM_JOB_ID}_ks_divergence"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"



### environment setup
export LANG=C
export LC_ALL=C
export PYTHONWARNINGS="ignore::FutureWarning"

eval "$(conda shell.bash hook)"



### define start time
START_TIME=$(date +%s)



### info
echo "Running on HYPERION - Date: $(date)"
echo "Tree file: ${TREE_FILE}"
echo "BUSCO sequences: ${BUSCO_SEQUENCES_DIR}"
echo "BUSCO runs: ${BUSCO_RUNS_DIR}"
echo "Results directory: ${RESULTS}"
echo "Outgroup: ${OUTGROUP}"



### sanity checks
if [ ! -f "${TREE_FILE}" ]; then
    echo "ERROR: Tree file not found: ${TREE_FILE}"
    exit 1
fi

if [ ! -d "${BUSCO_SEQUENCES_DIR}" ]; then
    echo "ERROR: BUSCO sequences directory not found: ${BUSCO_SEQUENCES_DIR}"
    exit 1
fi

if [ ! -d "${BUSCO_RUNS_DIR}" ]; then
    echo "ERROR: BUSCO runs directory not found: ${BUSCO_RUNS_DIR}"
    exit 1
fi

BUSCO_GENES_COUNT=$(ls -1 "${BUSCO_SEQUENCES_DIR}"/*.faa 2>/dev/null | wc -l)
echo "Found ${BUSCO_GENES_COUNT} BUSCO genes for analysis"



### ===========================================================================
### STEP 1: Parse tree and identify nodes for Ks analysis
### ===========================================================================
echo ""
echo "### STEP 1: Parsing phylogenetic tree ###"
echo ""

conda activate HTT_python_ete3

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

python3 << EOF
import sys
from pathlib import Path
from ete3 import Tree
import json

tree_file = Path("${TREE_FILE}")
results_dir = Path("${RESULTS}")
outgroup = "${OUTGROUP}"

tree = Tree(str(tree_file), format=1)

if outgroup in [leaf.name for leaf in tree.get_leaves()]:
    outgroup_node = tree.search_nodes(name=outgroup)[0]
    tree.set_outgroup(outgroup_node)
    print("Tree rooted with outgroup: " + outgroup)
else:
    print("Warning: Outgroup " + outgroup + " not found, assuming tree is already rooted")

nodes_info = []
node_counter = 0

for node in tree.traverse():
    if not node.is_leaf() and len(node.get_leaves()) > 1:
        node_counter += 1
        node_name = "node_" + str(node_counter)
        
        daughter_clades = []
        for child in node.get_children():
            clade_leaves = [leaf.name for leaf in child.get_leaves()]
            daughter_clades.append(clade_leaves)
        
        if len(daughter_clades) == 2 and any(outgroup in clade for clade in daughter_clades):
            outgroup_clade = [clade for clade in daughter_clades if outgroup in clade][0]
            if len(outgroup_clade) == 1:
                print("Skipping outgroup split: " + node_name)
                continue
        
        nodes_info.append({
            'node_name': node_name,
            'node_support': getattr(node, 'support', 0.0),
            'daughter_clades': daughter_clades,
            'all_leaves': [leaf.name for leaf in node.get_leaves()]
        })
        
        print("Node " + node_name + ": " + str(len(daughter_clades)) + " daughter clades, " + str(len(node.get_leaves())) + " total leaves")

print("Found " + str(len(nodes_info)) + " internal nodes for Ks analysis")

nodes_file = results_dir / "node_analyses" / "nodes_info.json"
with open(nodes_file, 'w') as f:
    json.dump(nodes_info, f, indent=2)

print("Node information saved to: " + str(nodes_file))
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed"
    exit 1
fi

echo "STEP 1 complete."
conda deactivate



### ===========================================================================
### STEP 2: For each node, select representative genes
### ===========================================================================
echo ""
echo "### STEP 2: Selecting representative genes ###"
echo ""

conda activate HTT_python_ete3

python3 << EOF
import json
import random
from pathlib import Path

results_dir = Path("${RESULTS}")
busco_dir = Path("${BUSCO_SEQUENCES_DIR}")

nodes_file = results_dir / "node_analyses" / "nodes_info.json"
with open(nodes_file) as f:
    nodes_info = json.load(f)

available_genes = [f.stem for f in busco_dir.glob("*.faa")]
print("Available BUSCO genes: " + str(len(available_genes)))

node_gene_pairs = []

for node_info in nodes_info:
    node_name = node_info['node_name']
    daughter_clades = node_info['daughter_clades']
    
    print("Processing " + node_name + " with " + str(len(daughter_clades)) + " daughter clades")
    
    for gene in available_genes:
        gene_file = busco_dir / (gene + ".faa")
        
        sequences = {}
        current_seq = ""
        current_name = ""
        
        with open(gene_file) as f:
            for line in f:
                if line.startswith(">"):
                    if current_name:
                        sequences[current_name] = current_seq
                    current_name = line[1:].strip()
                    current_seq = ""
                else:
                    current_seq += line.strip()
            if current_name:
                sequences[current_name] = current_seq
        
        representatives = []
        
        for i, clade in enumerate(daughter_clades):
            clade_seqs = {name: seq for name, seq in sequences.items() if name in clade}
            
            if not clade_seqs:
                break
            
            best_name = max(clade_seqs.keys(), key=lambda x: (len(clade_seqs[x]), random.random()))
            representatives.append((best_name, clade_seqs[best_name]))
        
        if len(representatives) == len(daughter_clades):
            node_gene_pairs.append({
                'node_name': node_name,
                'gene': gene,
                'representatives': representatives
            })

print("Total node-gene pairs for analysis: " + str(len(node_gene_pairs)))

pairs_file = results_dir / "node_analyses" / "node_gene_pairs.json"
with open(pairs_file, 'w') as f:
    json.dump(node_gene_pairs, f, indent=2)

print("Node-gene pairs saved to: " + str(pairs_file))
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 failed"
    exit 1
fi

echo "STEP 2 complete."
conda deactivate



### ===========================================================================
### STEP 3: Create protein alignments with MAFFT
### ===========================================================================
echo ""
echo "### STEP 3: Creating protein alignments with MAFFT ###"
echo ""

conda activate HTT_mafft

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

python3 << EOF
import json
import subprocess
from pathlib import Path

results_dir = Path("${RESULTS}")
pairs_file = results_dir / "node_analyses" / "node_gene_pairs.json"
alignments_dir = results_dir / "alignments"

with open(pairs_file) as f:
    node_gene_pairs = json.load(f)

print("Creating alignments for " + str(len(node_gene_pairs)) + " node-gene pairs")

success_count = 0
for i, pair in enumerate(node_gene_pairs):
    node_name = pair['node_name']
    gene = pair['gene']
    representatives = pair['representatives']
    
    if (i + 1) % 500 == 0:
        print("Progress: " + str(i + 1) + "/" + str(len(node_gene_pairs)))
    
    input_file = alignments_dir / (node_name + "_" + gene + "_input.faa")
    with open(input_file, 'w') as f:
        for name, seq in representatives:
            f.write(">" + name + "\n" + seq + "\n")
    
    output_file = alignments_dir / (node_name + "_" + gene + "_aligned.faa")
    
    try:
        with open(output_file, 'w') as out:
            subprocess.run([
                'mafft', '--auto', '--quiet',
                str(input_file)
            ], stdout=out, stderr=subprocess.DEVNULL, check=True)
        
        input_file.unlink()
        success_count += 1
        
    except subprocess.CalledProcessError:
        if input_file.exists():
            input_file.unlink()
        continue

print("Protein alignments complete: " + str(success_count) + " successful")
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Step 3 failed"
    exit 1
fi

echo "STEP 3 complete."
conda deactivate



### ===========================================================================
### STEP 4: Get nucleotide sequences from BUSCO and create codon alignments
### ===========================================================================
echo ""
echo "### STEP 4: Creating codon alignments with PAL2NAL ###"
echo ""

conda activate HTT_pal2nal

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
which pal2nal.pl

python3 << EOF
import json
import subprocess
from pathlib import Path

results_dir = Path("${RESULTS}")
pairs_file = results_dir / "node_analyses" / "node_gene_pairs.json"
alignments_dir = results_dir / "alignments"
codon_dir = results_dir / "codon_alignments"
busco_runs_dir = Path("${BUSCO_RUNS_DIR}")

codon_dir.mkdir(exist_ok=True)

with open(pairs_file) as f:
    node_gene_pairs = json.load(f)

print("Processing " + str(len(node_gene_pairs)) + " node-gene pairs for codon alignments")

# Build mapping: sample_name -> path to single_copy_busco_sequences (dothideomycetes)
sample_busco_dirs = {}
for busco_dir in busco_runs_dir.glob("*_dothideomycetes_odb10"):
    sample_name = busco_dir.name.replace("_dothideomycetes_odb10", "")
    seq_dir = busco_dir / "run_dothideomycetes_odb10" / "busco_sequences" / "single_copy_busco_sequences"
    if seq_dir.exists():
        sample_busco_dirs[sample_name] = seq_dir

print("Found " + str(len(sample_busco_dirs)) + " samples with BUSCO nucleotide sequences")

def read_fasta_first_seq(filepath):
    """Read first sequence from fasta file"""
    seq = ""
    with open(filepath) as f:
        for line in f:
            if line.startswith(">"):
                if seq:
                    return seq
            else:
                seq += line.strip()
    return seq

successful_pairs = 0
failed_pairs = 0

for i, pair in enumerate(node_gene_pairs):
    node_name = pair['node_name']
    gene = pair['gene']
    representatives = pair['representatives']
    
    if (i + 1) % 500 == 0:
        print("Progress: " + str(i + 1) + "/" + str(len(node_gene_pairs)) + " (success: " + str(successful_pairs) + ")")
    
    protein_alignment_file = alignments_dir / (node_name + "_" + gene + "_aligned.faa")
    if not protein_alignment_file.exists():
        failed_pairs += 1
        continue
    
    # Get nucleotide sequences
    nuc_sequences = {}
    all_found = True
    
    for rep_name, rep_seq in representatives:
        if rep_name not in sample_busco_dirs:
            all_found = False
            break
        
        nuc_file = sample_busco_dirs[rep_name] / (gene + ".fna")
        if not nuc_file.exists():
            all_found = False
            break
        
        nuc_seq = read_fasta_first_seq(nuc_file)
        if nuc_seq:
            nuc_sequences[rep_name] = nuc_seq
        else:
            all_found = False
            break
    
    if not all_found or len(nuc_sequences) != len(representatives):
        failed_pairs += 1
        continue
    
    # Write nucleotide sequences with headers matching protein alignment
    nuc_temp_file = codon_dir / (node_name + "_" + gene + "_nuc.fasta")
    with open(nuc_temp_file, 'w') as f:
        for name, seq in nuc_sequences.items():
            f.write(">" + name + "\n" + seq + "\n")
    
    codon_alignment_file = codon_dir / (node_name + "_" + gene + "_codon_aligned.fasta")
    
    try:
        with open(codon_alignment_file, 'w') as out:
            subprocess.run([
                'pal2nal.pl',
                str(protein_alignment_file),
                str(nuc_temp_file),
                '-output', 'fasta',
                '-nogap'
            ], stdout=out, stderr=subprocess.DEVNULL, check=True)
        
        if codon_alignment_file.stat().st_size > 0:
            successful_pairs += 1
        else:
            codon_alignment_file.unlink()
            failed_pairs += 1
        
    except subprocess.CalledProcessError:
        failed_pairs += 1
        if codon_alignment_file.exists():
            codon_alignment_file.unlink()
    
    if nuc_temp_file.exists():
        nuc_temp_file.unlink()

print("")
print("Successfully created " + str(successful_pairs) + " codon alignments")
print("Failed: " + str(failed_pairs))
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Step 4 failed"
    exit 1
fi

echo "STEP 4 complete."
conda deactivate



### ===========================================================================
### STEP 5: Calculate Ks values with R seqinr
### ===========================================================================
echo ""
echo "### STEP 5: Calculating Ks values ###"
echo ""

conda activate HTT_r

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

Rscript --vanilla - << EOF
library(seqinr)
library(jsonlite)

results_dir <- "${RESULTS}"
codon_dir <- file.path(results_dir, "codon_alignments")
ks_dir <- file.path(results_dir, "ks_calculations")
dir.create(ks_dir, showWarnings = FALSE)

codon_files <- list.files(codon_dir, pattern = "_codon_aligned.fasta\$", full.names = TRUE)
cat("Found", length(codon_files), "codon alignments\n")

if (length(codon_files) == 0) {
    cat("ERROR: No codon alignments found\n")
    quit(status = 1)
}

ks_results <- list()

for (i in seq_along(codon_files)) {
    codon_file <- codon_files[i]
    base_name <- basename(codon_file)
    pair_name <- sub("_codon_aligned.fasta\$", "", base_name)
    
    if (i %% 500 == 0) {
        cat("Progress:", i, "/", length(codon_files), "\n")
    }
    
    tryCatch({
        alignment <- read.alignment(codon_file, format = "fasta")
        
        if (length(alignment\$seq) == 2) {
            ks_result <- kaks(alignment)
            
            # Coerce to plain numeric (avoid 'dist' objects in JSON)
            ks_value <- as.numeric(ks_result\$ks)[1]
            ka_value <- as.numeric(ks_result\$ka)[1]
            
            node_gene <- strsplit(pair_name, "_", fixed = TRUE)[[1]]
            node_name <- paste(node_gene[1:2], collapse = "_")
            gene_name <- paste(node_gene[3:length(node_gene)], collapse = "_")
            
            ks_results[[pair_name]] <- list(
                node  = node_name,
                gene  = gene_name,
                ks    = ks_value,
                ka    = ka_value,
                ka_ks = ifelse(is.finite(ks_value) && ks_value > 0,
                               ka_value / ks_value,
                               NA_real_)
            )
        }
        
    }, error = function(e) {
        # silently skip errors
    })
}

cat("Calculated Ks for", length(ks_results), "gene pairs\n")

if (length(ks_results) == 0) {
    cat("ERROR: No Ks values calculated\n")
    quit(status = 1)
}

ks_file <- file.path(ks_dir, "ks_values.json")
write_json(ks_results, ks_file, pretty = TRUE)

# Calculate node-level statistics
node_stats <- list()

for (result in ks_results) {
    node   <- result\$node
    ks_val <- result\$ks
    
    if (is.finite(ks_val) && ks_val > 0) {
        if (is.null(node_stats[[node]])) {
            node_stats[[node]] <- c()
        }
        node_stats[[node]] <- c(node_stats[[node]], ks_val)
    }
}

node_summary <- list()

for (node in names(node_stats)) {
    ks_values <- node_stats[[node]]
    
    if (length(ks_values) >= 5) {
        node_summary[[node]] <- list(
            n_genes        = length(ks_values),
            ks_mean        = mean(ks_values),
            ks_median      = median(ks_values),
            ks_0.5_quantile = as.numeric(quantile(ks_values, 0.005)),
            ks_1_quantile   = as.numeric(quantile(ks_values, 0.01)),
            ks_5_quantile   = as.numeric(quantile(ks_values, 0.05)),
            ks_min         = min(ks_values),
            ks_max         = max(ks_values)
        )
        
        cat("Node", node, "- Genes:", length(ks_values), 
            "- Ks 0.5% quantile:", round(as.numeric(quantile(ks_values, 0.005)), 4), "\n")
    } else {
        cat("Node", node, "- Insufficient data (", length(ks_values), "genes)\n")
    }
}

summary_file <- file.path(ks_dir, "node_ks_summary.json")
write_json(node_summary, summary_file, pretty = TRUE)

if (length(node_summary) > 0) {
    summary_df <- do.call(rbind, lapply(names(node_summary), function(node) {
        stats <- node_summary[[node]]
        data.frame(
            node           = node,
            n_genes        = stats\$n_genes,
            ks_min         = stats\$ks_min,
            ks_0.5_quantile = stats\$ks_0.5_quantile,
            ks_1_quantile   = stats\$ks_1_quantile,
            ks_5_quantile   = stats\$ks_5_quantile,
            ks_median      = stats\$ks_median,
            ks_mean        = stats\$ks_mean,
            ks_max         = stats\$ks_max,
            stringsAsFactors = FALSE
        )
    }))
    
    summary_df\$pass_0.01 <- summary_df\$ks_0.5_quantile > 0.01
    summary_df\$pass_0.005 <- summary_df\$ks_0.5_quantile > 0.005
    summary_df\$pass_0.001 <- summary_df\$ks_0.5_quantile > 0.001
    
    write.table(summary_df, file.path(ks_dir, "node_ks_summary.tsv"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("\nNode summary:\n")
    print(summary_df)
}

# Save raw Ks values
ks_raw_df <- do.call(rbind, lapply(names(ks_results), function(pair) {
    data.frame(
        node = ks_results[[pair]]\$node,
        gene = ks_results[[pair]]\$gene,
        ks   = ks_results[[pair]]\$ks,
        ka   = ks_results[[pair]]\$ka,
        stringsAsFactors = FALSE
    )
}))
write.table(ks_raw_df, file.path(ks_dir, "ks_values_by_node.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nKs calculation complete\n")
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Step 5 failed"
    exit 1
fi

echo "STEP 5 complete."
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
echo "  Node information:     ${RESULTS}/node_analyses/nodes_info.json"
echo "  Node-gene pairs:      ${RESULTS}/node_analyses/node_gene_pairs.json"
echo "  Protein alignments:   ${RESULTS}/alignments/"
echo "  Codon alignments:     ${RESULTS}/codon_alignments/"
echo "  Ks values (JSON):     ${RESULTS}/ks_calculations/ks_values.json"
echo "  Ks by node (TSV):     ${RESULTS}/ks_calculations/ks_values_by_node.tsv"
echo "  Node Ks summary:      ${RESULTS}/ks_calculations/node_ks_summary.tsv"
echo ""
echo "Transfer to local machine for visualization:"
echo "ks_values_by_node.tsv & node_ks_summary.tsv & nodes_info.json"
echo ""