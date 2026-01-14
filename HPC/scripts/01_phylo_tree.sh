#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
GENOMES_DIR=${ROOT}/data/reference_genomes
RESULTS=${ROOT}/results/01_results_phylo_tree
LOGS=${ROOT}/logs/01_logs_phylo_tree
BUSCO_DB_PATH=/path/to/busco_downloads

### define outgroup for tree rooting
# IMPORTANT:
# The outgroup name must EXACTLY match the genome identifier (i.e. file prefix or assembly name) as provided in the downloaded genome dataset.
# Users MUST update this value to reflect the actual outgroup present in their local genome directory.
OUTGROUP="PsefloM405_1_AssemblyScaffolds_2024-02-18"



### create results/log folder(s)
mkdir -p "${RESULTS}"/{busco_runs,busco_summaries,sequences,alignments}
mkdir -p "${LOGS}"



### setup logging
JOBTAG="${SLURM_JOB_ID}_phylo_tree"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"



### environment setup
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export LANGUAGE=C
export PYTHONWARNINGS="ignore::FutureWarning"

eval "$(conda shell.bash hook)"



### define start time
START_TIME=$(date +%s)



### info
echo "Date: $(date)"
echo "Genomes directory: ${GENOMES_DIR}"
echo "Results directory: ${RESULTS}"
echo "BUSCO database path: ${BUSCO_DB_PATH}"
echo "Outgroup: ${OUTGROUP}"



### sanity checks
if [ ! -d "${GENOMES_DIR}" ]; then
    echo "ERROR: Genomes directory not found: ${GENOMES_DIR}"
    exit 1
fi

GENOME_COUNT=$(ls -1 "${GENOMES_DIR}"/*.fasta 2>/dev/null | wc -l)
if [ "${GENOME_COUNT}" -eq 0 ]; then
    echo "ERROR: No .fasta files found in ${GENOMES_DIR}"
    exit 1
fi
echo "Found ${GENOME_COUNT} genome assemblies"



### decompress any .gz files in genomes directory
echo "Checking for compressed genome FASTA files..."

shopt -s nullglob
for gz in "${GENOMES_DIR}"/*.gz; do
    echo "Decompressing: $(basename "${gz}")"
    gunzip -k "${gz}"
done
shopt -u nullglob



### define BUSCO lineages (taxonomically relevant for Cenococcum geophilum)
### Cenococcum: Dothideomycetes
LINEAGES=(
    "fungi_odb10"            # Kingdom
    "ascomycota_odb10"       # Phylum
    "dothideomycetes_odb10"  # Class
)



### ===========================================================================
### STEP 1: Run BUSCO for each genome with all lineages
### ===========================================================================
echo ""
echo "### STEP 1: Running BUSCO analysis (${#LINEAGES[@]} lineages per genome) ###"
echo ""

conda deactivate
conda activate HTT_busco

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "BUSCO version: $(busco --version)"

for genome in "${GENOMES_DIR}"/*.fasta; do
    GENOME_BASENAME=$(basename "${genome}")
    GENOME_NAME="${GENOME_BASENAME%%.*}"
    
    for lineage in "${LINEAGES[@]}"; do
        echo ""
        echo "Processing: ${GENOME_NAME} with lineage ${lineage}"
        
        OUTDIR="${RESULTS}/busco_runs/${GENOME_NAME}_${lineage}"
        
        if [ -d "${OUTDIR}" ]; then
            echo "  Skipping - output exists: ${OUTDIR}"
            continue
        fi
        
        busco \
            -i "${genome}" \
            -o "${GENOME_NAME}_${lineage}" \
            -l "${lineage}" \
            -m genome \
            -c "${SLURM_CPUS_PER_TASK}" \
            --out_path "${RESULTS}/busco_runs" \
            --download_path "${BUSCO_DB_PATH}" \
            --offline
        
        echo "  BUSCO complete for ${GENOME_NAME} with ${lineage}"
    done
done

echo ""
echo "STEP 1 complete."

conda deactivate



### ===========================================================================
### STEP 2: Generate BUSCO summary table
### ===========================================================================
echo ""
echo "### STEP 2: Generating BUSCO summary table ###"
echo ""

conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"

python3 << 'EOF'
import os
from pathlib import Path
import re

busco_dir = Path(os.environ.get('RESULTS', '.') + '/busco_runs')
summary_file = Path(os.environ.get('RESULTS', '.') + '/busco_summaries/busco_completeness_table.tsv')
summary_file.parent.mkdir(exist_ok=True)

lineages = [
    "fungi_odb10",
    "ascomycota_odb10",
    "dothideomycetes_odb10"
]

# find all samples
samples = set()
for d in busco_dir.iterdir():
    if d.is_dir():
        for lin in lineages:
            if d.name.endswith(f"_{lin}"):
                sample = d.name.replace(f"_{lin}", "")
                samples.add(sample)
                break

samples = sorted(samples)
print(f"Found {len(samples)} samples")

# parse BUSCO short summaries
results = {s: {} for s in samples}

for sample in samples:
    for lineage in lineages:
        summary_path = busco_dir / f"{sample}_{lineage}" / f"short_summary.specific.{lineage}.{sample}_{lineage}.txt"
        
        if not summary_path.exists():
            results[sample][lineage] = "NA"
            continue
        
        with open(summary_path) as f:
            content = f.read()
        
        # parse completeness percentage
        match = re.search(r"C:(\d+\.\d+)%", content)
        if match:
            results[sample][lineage] = match.group(1)
        else:
            results[sample][lineage] = "NA"

# write summary table
with open(summary_file, "w") as out:
    header = ["Sample"] + lineages
    out.write("\t".join(header) + "\n")
    
    for sample in samples:
        row = [sample] + [results[sample].get(lin, "NA") for lin in lineages]
        out.write("\t".join(row) + "\n")

print(f"Summary table written to: {summary_file}")
EOF

conda deactivate



### ===========================================================================
### STEP 3: Extract single-copy BUSCOs for phylogeny (dothideomycetes_odb10 only)
### ===========================================================================
echo ""
echo "### STEP 3: Extracting single-copy orthologs (dothideomycetes_odb10) ###"
echo ""

conda activate HTT_python

PHYLO_LINEAGE="dothideomycetes_odb10"
MISSING_THRESHOLD=0.10  # max 10% missing genomes per gene

export PHYLO_LINEAGE MISSING_THRESHOLD

python3 << 'EOF'
import os
from pathlib import Path
from collections import defaultdict

results_path = os.environ.get('RESULTS', '.')
busco_dir = Path(results_path + '/busco_runs')
seq_dir = Path(results_path + '/sequences')
seq_dir.mkdir(exist_ok=True)

lineage = os.environ.get('PHYLO_LINEAGE', 'dothideomycetes_odb10')
missing_threshold = float(os.environ.get('MISSING_THRESHOLD', '0.10'))

# find all samples with dothideomycetes_odb10 results
samples = []
for d in busco_dir.iterdir():
    if d.is_dir() and d.name.endswith(f"_{lineage}"):
        sample_name = d.name.replace(f"_{lineage}", "")
        samples.append(sample_name)

n_samples = len(samples)
min_present = int(n_samples * (1 - missing_threshold))
print(f"Found {n_samples} samples with {lineage} BUSCO results")
print(f"Requiring presence in at least {min_present} samples (>={100*(1-missing_threshold):.0f}%)")

# collect sequences per BUSCO
busco_seqs = defaultdict(dict)
total_buscos = set()

for sample in samples:
    sc_dir = busco_dir / f"{sample}_{lineage}" / f"run_{lineage}" / "busco_sequences" / "single_copy_busco_sequences"
    if not sc_dir.exists():
        print(f"  Warning: No single-copy sequences found for {sample}")
        continue
    
    for faa in sc_dir.glob("*.faa"):
        busco_id = faa.stem
        total_buscos.add(busco_id)
        with open(faa) as f:
            seq = "".join(line.strip() for line in f if not line.startswith(">"))
            busco_seqs[busco_id][sample] = seq

print(f"Total BUSCOs in {lineage}: {len(total_buscos)}")

# filter genes by presence threshold
selected = []
for busco_id, seqs in busco_seqs.items():
    if len(seqs) >= min_present:
        selected.append(busco_id)
        with open(seq_dir / f"{busco_id}.faa", "w") as out:
            for sample, seq in sorted(seqs.items()):
                out.write(f">{sample}\n{seq}\n")

print(f"Selected {len(selected)} / {len(total_buscos)} BUSCOs present in >= {min_present}/{n_samples} samples")
print(f"Sequence files written to: {seq_dir}")

# write gene list
with open(Path(results_path) / "selected_buscos.txt", "w") as f:
    for busco_id in sorted(selected):
        f.write(f"{busco_id}\n")
EOF

echo ""
echo "STEP 3 complete."

conda deactivate



### ===========================================================================
### STEP 4: Align sequences with FAMSA
### ===========================================================================
echo ""
echo "### STEP 4: Aligning sequences with FAMSA ###"
echo ""

conda activate HTT_famsa

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "FAMSA version: $(famsa 2>&1 | grep -i version | head -1)"

ALN_DIR="${RESULTS}/alignments"
mkdir -p "${ALN_DIR}"

for faa in "${RESULTS}"/sequences/*.faa; do
    BUSCO_ID=$(basename "${faa%.faa}")
    echo "Aligning: ${BUSCO_ID}"
    
    famsa -t "${SLURM_CPUS_PER_TASK}" "${faa}" "${ALN_DIR}/${BUSCO_ID}.aln"
done

echo ""
echo "STEP 4 complete."

conda deactivate



### ===========================================================================
### STEP 5: Trim alignments with trimAl
### ===========================================================================
echo ""
echo "### STEP 5: Trimming alignments with trimAl ###"
echo ""

conda activate HTT_trimal

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "trimAl version: $(trimal --version 2>&1 | head -1)"

for aln in "${ALN_DIR}"/*.aln; do
    # skip already trimmed files
    [[ "${aln}" == *".trimmed.aln" ]] && continue
    
    BUSCO_ID=$(basename "${aln%.aln}")
    echo "Trimming: ${BUSCO_ID}"
    
    trimal \
        -in "${aln}" \
        -out "${ALN_DIR}/${BUSCO_ID}.trimmed.aln" \
        -automated1
done

echo ""
echo "STEP 5 complete."

conda deactivate



### ===========================================================================
### STEP 6: Concatenate alignments
### ===========================================================================
echo ""
echo "### STEP 6: Concatenating alignments ###"
echo ""

conda activate HTT_python

export ALN_DIR

python3 << 'EOF'
import os
from pathlib import Path
from collections import defaultdict

aln_dir = Path(os.environ.get('ALN_DIR', '.'))
results_path = os.environ.get('RESULTS', '.')
out_file = Path(results_path) / "concatenated.faa"
partition_file = Path(results_path) / "partitions.txt"

samples = set()
gene_alns = {}

# parse all trimmed alignments
for aln in sorted(aln_dir.glob("*.trimmed.aln")):
    gene = aln.stem.replace(".trimmed", "")
    seqs = {}
    current = None
    
    with open(aln) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                samples.add(current)
                seqs[current] = []
            else:
                seqs[current].append(line.strip())
    
    gene_alns[gene] = {k: "".join(v) for k, v in seqs.items()}

samples = sorted(samples)
genes = sorted(gene_alns.keys())

print(f"Concatenating {len(genes)} genes for {len(samples)} samples")

# write concatenated alignment
with open(out_file, "w") as out:
    for sample in samples:
        out.write(f">{sample}\n")
        concat = []
        for gene in genes:
            if sample in gene_alns[gene]:
                concat.append(gene_alns[gene][sample])
            else:
                # fill with gaps for missing data
                length = len(next(iter(gene_alns[gene].values())))
                concat.append("-" * length)
        out.write("".join(concat) + "\n")

# write partition file (for gene-wise model testing)
pos = 1
with open(partition_file, "w") as pf:
    for gene in genes:
        length = len(next(iter(gene_alns[gene].values())))
        pf.write(f"AUTO, {gene} = {pos}-{pos + length - 1}\n")
        pos += length

print(f"Concatenated alignment: {out_file}")
print(f"Partition file: {partition_file}")
print(f"Total alignment length: {pos - 1} aa")
EOF

echo ""
echo "STEP 6 complete."

conda deactivate



### ===========================================================================
### STEP 7: Build phylogenetic tree with IQ-TREE
### ===========================================================================
echo ""
echo "### STEP 7: Building phylogenetic tree with IQ-TREE ###"
echo ""

conda activate HTT_iqtree

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "IQ-TREE version: $(iqtree2 --version 2>&1 | head -1)"

cd "${RESULTS}"

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
if (!outgroup %in% tree\$tip.label) {
    cat("ERROR: Outgroup not found in tree\n")
    cat("Looking for:", outgroup, "\n")
    cat("Available tips:\n")
    cat(paste(tree\$tip.label, collapse = "\n"), "\n")
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
echo "  BUSCO results:        ${RESULTS}/busco_runs/"
echo "  BUSCO summary table:  ${RESULTS}/busco_summaries/busco_completeness_table.tsv"
echo "  Selected BUSCOs:      ${RESULTS}/selected_buscos.txt"
echo "  Sequences:            ${RESULTS}/sequences/"
echo "  Alignments:           ${RESULTS}/alignments/"
echo "  Concatenated:         ${RESULTS}/concatenated.faa"
echo "  Partitions:           ${RESULTS}/partitions.txt"
echo "  Tree file (unrooted): ${RESULTS}/phylogeny.treefile"
echo "  Tree file (rooted):   ${RESULTS}/phylogeny_rooted.treefile"
echo "  IQ-TREE log:          ${RESULTS}/phylogeny.iqtree"
echo ""