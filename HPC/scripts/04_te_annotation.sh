#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
GENOMES_DIR=${ROOT}/data/reference_genomes
RESULTS=${ROOT}/results/04_results_te_annotation
LOGS=${ROOT}/logs/04_logs_te_annotation

### TE filtering parameters
# IMPORTANT:
# These parameters control which TEs are retained for downstream HTT analysis.
# Users may adjust based on their specific research questions.
MIN_TE_LENGTH=300  # minimum TE length in bp



### create results/log folder(s)
mkdir -p "${RESULTS}"/{earlgrey_runs,te_sequences,te_filtered,summary}
mkdir -p "${LOGS}"



### setup logging
# supports both array jobs and single jobs
if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
    JOBTAG="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_te_annotation"
else
    JOBTAG="${SLURM_JOB_ID}_te_annotation"
fi
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
echo "Minimum TE length: ${MIN_TE_LENGTH} bp"



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



### determine run mode (array job vs single job)
if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
    # array job mode: process single genome
    RUN_MODE="array"
    GENOMES=(${GENOMES_DIR}/*.fasta)
    
    if [ ${SLURM_ARRAY_TASK_ID} -ge ${GENOME_COUNT} ]; then
        echo "Task ID ${SLURM_ARRAY_TASK_ID} exceeds genome count ${GENOME_COUNT}, exiting"
        exit 0
    fi
    
    GENOME_LIST=("${GENOMES[${SLURM_ARRAY_TASK_ID}]}")
    echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
    echo "Processing genome: $(basename ${GENOME_LIST[0]})"
else
    # single job mode: process all genomes sequentially
    RUN_MODE="single"
    GENOME_LIST=(${GENOMES_DIR}/*.fasta)
    echo "Single job mode: processing all ${GENOME_COUNT} genomes"
fi



### ===========================================================================
### STEP 1: Run EarlGrey TE annotation
### ===========================================================================
echo ""
echo "### STEP 1: Running EarlGrey TE annotation ###"
echo ""

conda deactivate
conda activate HTT_earlgrey

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "EarlGrey version: $(earlGrey -v 2>&1 | head -1)"

for GENOME in "${GENOME_LIST[@]}"; do
    GENOME_BASENAME=$(basename "${GENOME}" .fasta)
    
    # extract strain ID (handles both Cenge and Pseflo patterns)
    STRAIN_ID=$(echo "${GENOME_BASENAME}" | grep -oP '^(Cenge|Pseflo)[A-Za-z0-9]+')
    
    if [ -z "${STRAIN_ID}" ]; then
        echo "Warning: Could not extract strain ID from ${GENOME_BASENAME}, using full name"
        STRAIN_ID="${GENOME_BASENAME}"
    fi
    
    echo ""
    echo "Processing: ${STRAIN_ID}"
    
    # check for existing EarlGrey output
    EARLGREY_DIR=$(ls -d "${RESULTS}/earlgrey_runs/${STRAIN_ID}_EarlGrey"* 2>/dev/null | head -1)
    
    if [ -n "${EARLGREY_DIR}" ] && [ -d "${EARLGREY_DIR}" ]; then
        echo "  Skipping - output exists: ${EARLGREY_DIR}"
    else
        earlGrey \
            -g "${GENOME}" \
            -s "${STRAIN_ID}" \
            -o "${RESULTS}/earlgrey_runs" \
            -t "${SLURM_CPUS_PER_TASK}"
        
        echo "  EarlGrey complete for ${STRAIN_ID}"
    fi
done

echo ""
echo "STEP 1 complete."

conda deactivate



### ===========================================================================
### STEP 2: Extract and filter TE sequences
### ===========================================================================
echo ""
echo "### STEP 2: Extracting and filtering TE sequences ###"
echo ""

conda activate HTT_bedtools

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "BEDtools version: $(bedtools --version)"

for GENOME in "${GENOME_LIST[@]}"; do
    GENOME_BASENAME=$(basename "${GENOME}" .fasta)
    
    # extract strain ID
    STRAIN_ID=$(echo "${GENOME_BASENAME}" | grep -oP '^(Cenge|Pseflo)[A-Za-z0-9]+')
    if [ -z "${STRAIN_ID}" ]; then
        STRAIN_ID="${GENOME_BASENAME}"
    fi
    
    echo ""
    echo "Processing: ${STRAIN_ID}"
    
    # find EarlGrey output directory
    EARLGREY_DIR=$(ls -d "${RESULTS}/earlgrey_runs/${STRAIN_ID}_EarlGrey"* 2>/dev/null | head -1)
    
    if [ -z "${EARLGREY_DIR}" ] || [ ! -d "${EARLGREY_DIR}" ]; then
        echo "  ERROR: No EarlGrey output found for ${STRAIN_ID}"
        continue
    fi
    
    # find GFF file
    GFF_FILE=$(find "${EARLGREY_DIR}" -name "*.filteredRepeats.gff" | head -1)
    
    if [ ! -f "${GFF_FILE}" ]; then
        echo "  ERROR: GFF file not found in ${EARLGREY_DIR}"
        continue
    fi
    
    echo "  Found GFF: $(basename ${GFF_FILE})"
    
    # output files
    FILTERED_GFF="${RESULTS}/te_filtered/${STRAIN_ID}_te_filtered.gff"
    TE_FASTA="${RESULTS}/te_sequences/${STRAIN_ID}_te.fasta"
    
    # filter TEs:
    # - Remove Satellite, Simple_repeat, Low_complexity, rRNA, tRNA, snRNA, scRNA
    # - Keep only TEs >= MIN_TE_LENGTH bp
    awk -F'\t' -v min_len="${MIN_TE_LENGTH}" '
    BEGIN {OFS="\t"}
    {
        # skip header lines
        if ($0 ~ /^#/) next
        
        # calculate length
        len = $5 - $4 + 1
        
        # skip if below minimum length
        if (len < min_len) next
        
        # get classification from column 3
        class = $3
        
        # skip non-TE elements
        if (class ~ /Satellite/) next
        if (class ~ /Simple_repeat/) next
        if (class ~ /Low_complexity/) next
        if (class ~ /rRNA/) next
        if (class ~ /tRNA/) next
        if (class ~ /snRNA/) next
        if (class ~ /scRNA/) next
        if (class ~ /-rich/) next
        
        # print passing TEs
        print $0
    }
    ' "${GFF_FILE}" > "${FILTERED_GFF}"
    
    # count filtered TEs
    TE_COUNT=$(wc -l < "${FILTERED_GFF}")
    echo "  Filtered TEs: ${TE_COUNT}"
    
    # create BED file for sequence extraction
    awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/{
        name=$3"::"$1":"$4"-"$5
        score=$6; if(score=="." || score=="") score=0
        print $1, $4-1, $5, name, score, $7
    }' "${FILTERED_GFF}" > "${FILTERED_GFF%.gff}.bed"
    
    # extract TE sequences
    bedtools getfasta \
        -fi "${GENOME}" \
        -bed "${FILTERED_GFF%.gff}.bed" \
        -fo "${TE_FASTA}" \
        -name \
        -s
    
    rm -f "${FILTERED_GFF%.gff}.bed"
    
    echo "  TE sequences: ${TE_FASTA}"
done

echo ""
echo "STEP 2 complete."

conda deactivate



### ===========================================================================
### STEP 3: Generate summary statistics (single job mode only, or final array task)
### ===========================================================================

# only run summary if single job mode or if SUMMARIZE=true is set
if [ "${RUN_MODE}" = "single" ] || [ "${SUMMARIZE}" = "true" ]; then

echo ""
echo "### STEP 3: Generating summary statistics ###"
echo ""

conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"

export RESULTS MIN_TE_LENGTH

python3 << 'EOF'
import os
from pathlib import Path

results_dir = Path(os.environ.get('RESULTS', '.'))
min_te_length = os.environ.get('MIN_TE_LENGTH', '300')
earlgrey_runs = results_dir / "earlgrey_runs"
te_filtered_dir = results_dir / "te_filtered"
te_sequences_dir = results_dir / "te_sequences"
summary_dir = results_dir / "summary"
summary_dir.mkdir(exist_ok=True)

# collect statistics per genome
stats = []
total_tes = 0
total_raw = 0
total_filtered_out = 0

for gff_file in sorted(te_filtered_dir.glob("*_te_filtered.gff")):
    genome_name = gff_file.stem.replace("_te_filtered", "")
    
    # find EarlGrey directory
    earlgrey_dirs = list(earlgrey_runs.glob(f"{genome_name}_EarlGrey*"))
    raw_count = 0
    
    if earlgrey_dirs:
        earlgrey_dir = earlgrey_dirs[0]
        raw_gff_list = list(earlgrey_dir.glob("*_summaryFiles/*.filteredRepeats.gff"))
        if raw_gff_list:
            raw_gff = raw_gff_list[0]
            with open(raw_gff) as f:
                for line in f:
                    if not line.startswith("#") and len(line.strip().split("\t")) >= 9:
                        raw_count += 1
    
    # count filtered TEs and classify
    te_classes = {}
    te_count = 0
    total_length = 0
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            te_count += 1
            length = int(fields[4]) - int(fields[3]) + 1
            total_length += length
            
            # TE class is in column 3
            te_class = fields[2].split("/")[0]
            te_classes[te_class] = te_classes.get(te_class, 0) + 1
    
    filtered_out = raw_count - te_count
    pct_removed = (filtered_out / raw_count * 100) if raw_count > 0 else 0
    
    stats.append({
        'genome': genome_name,
        'raw_count': raw_count,
        'te_count': te_count,
        'filtered_out': filtered_out,
        'pct_removed': pct_removed,
        'total_bp': total_length,
        'classes': te_classes
    })
    total_tes += te_count
    total_raw += raw_count
    total_filtered_out += filtered_out
    
    print(f"{genome_name}: {raw_count} -> {te_count} TEs (removed {filtered_out}, {pct_removed:.1f}%), {total_length:,} bp")

# write summary table
summary_file = summary_dir / "te_annotation_summary.tsv"
with open(summary_file, "w") as out:
    out.write("genome\tte_count\ttotal_bp\n")
    for s in stats:
        out.write(f"{s['genome']}\t{s['te_count']}\t{s['total_bp']}\n")

print(f"\nTotal TEs across all genomes: {total_tes:,}")
print(f"Summary written to: {summary_file}")

# write filtering summary
filter_file = summary_dir / "te_filtering_summary.tsv"
with open(filter_file, "w") as out:
    out.write(f"# Filtering parameters: min_length >= {min_te_length} bp\n")
    out.write("genome\traw_count\tfiltered_count\tremoved_count\tpct_removed\n")
    for s in stats:
        out.write(f"{s['genome']}\t{s['raw_count']}\t{s['te_count']}\t{s['filtered_out']}\t{s['pct_removed']:.1f}\n")
    # totals row
    total_pct = (total_filtered_out / total_raw * 100) if total_raw > 0 else 0
    out.write(f"TOTAL\t{total_raw}\t{total_tes}\t{total_filtered_out}\t{total_pct:.1f}\n")

print(f"Filtering summary written to: {filter_file}")
print(f"\nOverall: {total_raw:,} -> {total_tes:,} TEs (removed {total_filtered_out:,}, {total_pct:.1f}%)")

# write class breakdown
class_file = summary_dir / "te_class_summary.tsv"
all_classes = set()
for s in stats:
    all_classes.update(s['classes'].keys())
all_classes = sorted(all_classes)

with open(class_file, "w") as out:
    out.write("genome\t" + "\t".join(all_classes) + "\n")
    for s in stats:
        row = [s['genome']] + [str(s['classes'].get(c, 0)) for c in all_classes]
        out.write("\t".join(row) + "\n")

print(f"Class breakdown written to: {class_file}")
EOF

echo ""
echo "STEP 3 complete."

conda deactivate

fi



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
echo "  EarlGrey runs:        ${RESULTS}/earlgrey_runs/"
echo "  Filtered TE GFFs:     ${RESULTS}/te_filtered/"
echo "  TE sequences:         ${RESULTS}/te_sequences/"
echo "  Summary statistics:   ${RESULTS}/summary/te_annotation_summary.tsv"
echo "  Filtering summary:    ${RESULTS}/summary/te_filtering_summary.tsv"
echo "  Class breakdown:      ${RESULTS}/summary/te_class_summary.tsv"
echo ""
echo "Usage notes:"
echo "  - For array jobs, run summary separately after all tasks complete:"
echo "    SUMMARIZE=true sbatch 04_te_annotation.sh"
echo "  - Or run as single job to process all genomes sequentially"
echo ""