#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
GENOMES_DIR=${ROOT}/data/reference_genomes
RESULTS=${ROOT}/results/04_results_te_annotation
LOGS=${ROOT}/logs/04_logs_te_annotation



### create results/log folder(s)
mkdir -p "${RESULTS}"/{earlgrey_runs,te_sequences,te_filtered,summary}
mkdir -p "${LOGS}"



### setup logging
JOBTAG="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_te_annotation"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"



### environment setup
export LANG=C
export LC_ALL=C
unset LANGUAGE
export PYTHONWARNINGS="ignore::FutureWarning"

eval "$(conda shell.bash hook)"



### define start time
START_TIME=$(date +%s)



### get genome for this array task
GENOMES=(${GENOMES_DIR}/*.fasta)
GENOME_COUNT=${#GENOMES[@]}

# check if this task ID is valid
if [ ${SLURM_ARRAY_TASK_ID} -ge ${GENOME_COUNT} ]; then
    echo "Task ID ${SLURM_ARRAY_TASK_ID} exceeds genome count ${GENOME_COUNT}, exiting"
    exit 0
fi

GENOME="${GENOMES[${SLURM_ARRAY_TASK_ID}]}"
GENOME_BASENAME=$(basename "${GENOME}" .fasta)

# Extract strain ID (e.g., Cenge1005 from Cenge1005_1_AssemblyScaffolds_2024-05-13)
STRAIN_ID=$(echo "${GENOME_BASENAME}" | grep -oP '^Cenge\d+')

if [ -z "${STRAIN_ID}" ]; then
    echo "ERROR: Could not extract strain ID from ${GENOME_BASENAME}"
    exit 1
fi



### info
echo "Running on HYPERION - Date: $(date)"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Genome file: ${GENOME_BASENAME}"
echo "Strain ID: ${STRAIN_ID}"
echo "Genomes directory: ${GENOMES_DIR}"
echo "Results directory: ${RESULTS}"
echo "Total genomes: ${GENOME_COUNT}"



### sanity checks
if [ ! -f "${GENOME}" ]; then
    echo "ERROR: Genome file not found: ${GENOME}"
    exit 1
fi



### ===========================================================================
### STEP 1: Run EarlGrey TE annotation (or skip if exists)
### ===========================================================================
echo ""
echo "### STEP 1: EarlGrey TE annotation for ${STRAIN_ID} ###"
echo ""

# Check for existing EarlGrey output
EARLGREY_DIR=$(ls -d "${RESULTS}/earlgrey_runs/${STRAIN_ID}_EarlGrey"* 2>/dev/null | head -1)

if [ -n "${EARLGREY_DIR}" ] && [ -d "${EARLGREY_DIR}" ]; then
    echo "Skipping EarlGrey - using existing output: ${EARLGREY_DIR}"
else
    # Run EarlGrey
    conda activate HTT_earlgrey
    echo "Running in environment: ${CONDA_DEFAULT_ENV}"
    
    earlGrey \
        -g "${GENOME}" \
        -s "${STRAIN_ID}" \
        -o "${RESULTS}/earlgrey_runs" \
        -t "${SLURM_CPUS_PER_TASK}"
    
    echo "EarlGrey complete for ${STRAIN_ID}"
    conda deactivate
    
    # Find the newly created directory
    EARLGREY_DIR=$(ls -d "${RESULTS}/earlgrey_runs/${STRAIN_ID}_EarlGrey"* 2>/dev/null | head -1)
fi

# Verify EarlGrey directory exists
if [ -z "${EARLGREY_DIR}" ] || [ ! -d "${EARLGREY_DIR}" ]; then
    echo "ERROR: No EarlGrey output found for ${STRAIN_ID}"
    echo "Expected pattern: ${RESULTS}/earlgrey_runs/${STRAIN_ID}_EarlGrey*"
    exit 1
fi



### ===========================================================================
### STEP 2: Extract TE sequences and apply filters
### ===========================================================================
echo ""
echo "### STEP 2: Extracting and filtering TE sequences for ${STRAIN_ID} ###"
echo ""

conda activate HTT_bedtools

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "BEDtools version: $(bedtools --version)"

# Find the GFF file from EarlGrey output
GFF_FILE=$(find "${EARLGREY_DIR}" -name "*.filteredRepeats.gff" | head -1)

if [ ! -f "${GFF_FILE}" ]; then
    echo "ERROR: GFF file not found in: ${EARLGREY_DIR}"
    echo "EarlGrey may have failed for this genome"
    exit 1
fi

echo "Found GFF: ${GFF_FILE}"

# Output files
FILTERED_GFF="${RESULTS}/te_filtered/${STRAIN_ID}_te_filtered.gff"
TE_FASTA="${RESULTS}/te_sequences/${STRAIN_ID}_te.fasta"

# Filter TEs according to paper criteria:
# - Remove Satellite
# - Keep only TEs >= 300bp
# - Classification is in column 3
awk -F'\t' '
BEGIN {OFS="\t"}
{
    # skip header lines
    if ($0 ~ /^#/) next
    
    # calculate length
    len = $5 - $4 + 1
    
    # skip if < 300bp
    if (len < 300) next
    
    # get classification from column 3
    class = $3
    
    # skip satellite
    if (class ~ /Satellite/) next
    
    # skip simple repeats, low complexity if present
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

# Count filtered TEs
TE_COUNT=$(wc -l < "${FILTERED_GFF}")
echo "Filtered TEs: ${TE_COUNT}"

# Extract TE sequences with bedtools
awk -F'\t' 'BEGIN{OFS="\t"} $0 !~ /^#/{
  name=$3"::"$1":"$4"-"$5
  score=$6; if(score=="." || score=="") score=0
  print $1, $4-1, $5, name, score, $7
}' "${FILTERED_GFF}" > "${FILTERED_GFF%.gff}.bed"

bedtools getfasta \
  -fi "${GENOME}" \
  -bed "${FILTERED_GFF%.gff}.bed" \
  -fo "${TE_FASTA}" \
  -name \
  -s

rm -f "${FILTERED_GFF%.gff}.bed"

echo "TE sequences extracted to: ${TE_FASTA}"

conda deactivate



### ===========================================================================
### FINISH
### ===========================================================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "=========================================="
echo "Array task ${SLURM_ARRAY_TASK_ID} complete: ${STRAIN_ID}"
echo "Total runtime: $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"
echo "=========================================="
echo ""
echo "Output files:"
echo "  EarlGrey run:         ${EARLGREY_DIR}"
echo "  Filtered TE GFF:      ${FILTERED_GFF}"
echo "  TE sequences:         ${TE_FASTA}"
echo ""