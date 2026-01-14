#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
TE_RESULTS=${ROOT}/results/04_results_te_annotation
RESULTS=${ROOT}/results/05_results_te_classification
LOGS=${ROOT}/logs/05_logs_te_classification

### database paths
# IMPORTANT:
# Users MUST update these paths to point to their local database installations.
# CDD: Download from NCBI CDD (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/)
# Pfam: Download TE-specific HMMs or full Pfam-A.hmm from Pfam
CDD_DB=${ROOT}/databases/cdd/Cdd
PFAM_DB=${ROOT}/databases/pfam/Pfam-A.TEonly.hmm
TE_DOMAINS_LOOKUP=${ROOT}/data/te_domains_lookup.tsv

### search parameters
EVALUE_THRESHOLD=0.001



### create results/log folder(s)
mkdir -p "${RESULTS}"/{translated,cdd_hits,pfam_hits,filtered_tes,header_maps,deduplicated,summary}
mkdir -p "${LOGS}"



### setup logging
if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
    JOBTAG="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_te_classification"
else
    JOBTAG="${SLURM_JOB_ID}_te_classification"
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
echo "TE results directory: ${TE_RESULTS}"
echo "Results directory: ${RESULTS}"
echo "CDD database: ${CDD_DB}"
echo "Pfam database: ${PFAM_DB}"
echo "TE domains lookup: ${TE_DOMAINS_LOOKUP}"
echo "E-value threshold: ${EVALUE_THRESHOLD}"



### sanity checks
if [ ! -d "${TE_RESULTS}/te_sequences" ]; then
    echo "ERROR: TE sequences directory not found: ${TE_RESULTS}/te_sequences"
    echo "Run 04_te_annotation.sh first"
    exit 1
fi

if [ ! -f "${CDD_DB}.rps" ] && [ ! -f "${CDD_DB}.pal" ]; then
    echo "ERROR: CDD database not found: ${CDD_DB}"
    exit 1
fi

if [ ! -f "${PFAM_DB}" ]; then
    echo "ERROR: Pfam database not found: ${PFAM_DB}"
    exit 1
fi

if [ ! -f "${TE_DOMAINS_LOOKUP}" ]; then
    echo "ERROR: TE domains lookup file not found: ${TE_DOMAINS_LOOKUP}"
    exit 1
fi

TE_FASTA_COUNT=$(ls -1 "${TE_RESULTS}/te_sequences/"*_te.fasta 2>/dev/null | wc -l)
if [ "${TE_FASTA_COUNT}" -eq 0 ]; then
    echo "ERROR: No TE fasta files found in ${TE_RESULTS}/te_sequences/"
    exit 1
fi
echo "Found ${TE_FASTA_COUNT} TE fasta files"



### determine run mode (array job vs single job)
mapfile -t TE_FASTAS < <(ls -1 "${TE_RESULTS}/te_sequences/"*_te.fasta | sort)

if [ -n "${SLURM_ARRAY_TASK_ID}" ]; then
    RUN_MODE="array"
    
    if [ "${SLURM_ARRAY_TASK_ID}" -ge "${TE_FASTA_COUNT}" ]; then
        echo "Task ID ${SLURM_ARRAY_TASK_ID} exceeds TE fasta count ${TE_FASTA_COUNT}, exiting"
        exit 0
    fi
    
    TE_FASTA_LIST=("${TE_FASTAS[${SLURM_ARRAY_TASK_ID}]}")
    echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
    echo "Processing: $(basename ${TE_FASTA_LIST[0]})"
else
    RUN_MODE="single"
    TE_FASTA_LIST=("${TE_FASTAS[@]}")
    echo "Single job mode: processing all ${TE_FASTA_COUNT} genomes"
fi



### ===========================================================================
### STEP 1: Create simplified headers for transeq
### ===========================================================================
echo ""
echo "### STEP 1: Creating simplified headers for transeq ###"
echo ""

# transeq truncates long headers, so we create a mapping
for TE_FASTA in "${TE_FASTA_LIST[@]}"; do
    STRAIN_ID=$(basename "${TE_FASTA}" _te.fasta)
    
    SIMPLIFIED_FASTA="${RESULTS}/translated/${STRAIN_ID}_te_simplified.fasta"
    HEADER_MAP="${RESULTS}/header_maps/${STRAIN_ID}_header_map.tsv"
    
    if [ -f "${HEADER_MAP}" ]; then
        echo "  ${STRAIN_ID}: header map exists, skipping"
        continue
    fi
    
    echo "  ${STRAIN_ID}: creating simplified headers"
    
    awk '
    BEGIN { n=0 }
    /^>/ {
        n++
        full = substr($0, 2)
        split(full, arr, " ")
        full = arr[1]
        short = "TE" n
        print short "\t" full > "'"${HEADER_MAP}"'"
        print ">" short
        next
    }
    { print }
    ' "${TE_FASTA}" > "${SIMPLIFIED_FASTA}"
done

echo ""
echo "STEP 1 complete."



### ===========================================================================
### STEP 2: Translate sequences with transeq (6-frame)
### ===========================================================================
echo ""
echo "### STEP 2: Translating sequences with transeq ###"
echo ""

conda deactivate
conda activate HTT_emboss

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "transeq version: $(transeq -version 2>&1 | head -1)"

for TE_FASTA in "${TE_FASTA_LIST[@]}"; do
    STRAIN_ID=$(basename "${TE_FASTA}" _te.fasta)
    
    SIMPLIFIED_FASTA="${RESULTS}/translated/${STRAIN_ID}_te_simplified.fasta"
    TRANSLATED="${RESULTS}/translated/${STRAIN_ID}_te_translated.faa"
    
    if [ -f "${TRANSLATED}" ]; then
        echo "  ${STRAIN_ID}: translated file exists, skipping"
        continue
    fi
    
    echo "  ${STRAIN_ID}: running transeq"
    
    transeq \
        -sequence "${SIMPLIFIED_FASTA}" \
        -outseq "${TRANSLATED}" \
        -frame 6 \
        -clean
done

echo ""
echo "STEP 2 complete."

conda deactivate



### ===========================================================================
### STEP 3: Search CDD with rpstblastn
### ===========================================================================
echo ""
echo "### STEP 3: Searching CDD with rpstblastn ###"
echo ""

conda activate HTT_blast

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "rpstblastn version: $(rpstblastn -version 2>&1 | head -1)"

for TE_FASTA in "${TE_FASTA_LIST[@]}"; do
    STRAIN_ID=$(basename "${TE_FASTA}" _te.fasta)
    
    CDD_OUT="${RESULTS}/cdd_hits/${STRAIN_ID}_cdd_hits.tsv"
    
    if [ -f "${CDD_OUT}" ]; then
        echo "  ${STRAIN_ID}: CDD hits exist, skipping"
        continue
    fi
    
    echo "  ${STRAIN_ID}: running rpstblastn"
    
    rpstblastn \
        -query "${TE_FASTA}" \
        -db "${CDD_DB}" \
        -out "${CDD_OUT}" \
        -evalue "${EVALUE_THRESHOLD}" \
        -max_target_seqs 1 \
        -num_threads "${SLURM_CPUS_PER_TASK}" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qframe"
done

echo ""
echo "STEP 3 complete."

conda deactivate



### ===========================================================================
### STEP 4: Search Pfam with hmmsearch
### ===========================================================================
echo ""
echo "### STEP 4: Searching Pfam with hmmsearch ###"
echo ""

conda activate HTT_hmmer

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "hmmsearch version: $(hmmsearch -h 2>&1 | grep "HMMER" | head -1)"

for TE_FASTA in "${TE_FASTA_LIST[@]}"; do
    STRAIN_ID=$(basename "${TE_FASTA}" _te.fasta)
    
    TRANSLATED="${RESULTS}/translated/${STRAIN_ID}_te_translated.faa"
    PFAM_OUT="${RESULTS}/pfam_hits/${STRAIN_ID}_pfam_hits.tbl"
    PFAM_DOM="${RESULTS}/pfam_hits/${STRAIN_ID}_pfam_domtbl.tbl"
    
    if [ -f "${PFAM_DOM}" ]; then
        echo "  ${STRAIN_ID}: Pfam hits exist, skipping"
        continue
    fi
    
    echo "  ${STRAIN_ID}: running hmmsearch"
    
    hmmsearch \
        --tblout "${PFAM_OUT}" \
        --domtblout "${PFAM_DOM}" \
        -E "${EVALUE_THRESHOLD}" \
        --cpu "${SLURM_CPUS_PER_TASK}" \
        "${PFAM_DB}" \
        "${TRANSLATED}"
done

echo ""
echo "STEP 4 complete."

conda deactivate



### ===========================================================================
### STEP 5: Filter TEs by domain hits
### ===========================================================================
echo ""
echo "### STEP 5: Filtering TEs by domain hits ###"
echo ""

conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"

for TE_FASTA in "${TE_FASTA_LIST[@]}"; do
    STRAIN_ID=$(basename "${TE_FASTA}" _te.fasta)
    
    CDD_OUT="${RESULTS}/cdd_hits/${STRAIN_ID}_cdd_hits.tsv"
    PFAM_DOM="${RESULTS}/pfam_hits/${STRAIN_ID}_pfam_domtbl.tbl"
    HEADER_MAP="${RESULTS}/header_maps/${STRAIN_ID}_header_map.tsv"
    
    echo "  ${STRAIN_ID}: filtering by domain hits"
    
    export TE_FASTA CDD_OUT PFAM_DOM RESULTS TE_DOMAINS_LOOKUP STRAIN_ID HEADER_MAP
    
    python3 << 'EOF'
import os
from pathlib import Path
from collections import defaultdict
import re

te_fasta = Path(os.environ.get('TE_FASTA', ''))
cdd_file = Path(os.environ.get('CDD_OUT', ''))
pfam_file = Path(os.environ.get('PFAM_DOM', ''))
results_dir = Path(os.environ.get('RESULTS', '.'))
lookup_file = Path(os.environ.get('TE_DOMAINS_LOOKUP', ''))
strain = os.environ.get('STRAIN_ID', '')
header_map_file = Path(os.environ.get('HEADER_MAP', ''))

filtered_dir = results_dir / "filtered_tes"
filtered_dir.mkdir(exist_ok=True)

# load header mapping (short -> full)
short_to_full = {}
with open(header_map_file) as f:
    for line in f:
        short, full = line.rstrip("\n").split("\t")
        short_to_full[short] = full

# read original fasta to get TE lengths
def read_fasta_lengths(fp: Path):
    lengths = {}
    cur = None
    L = 0
    for line in fp.open():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur is not None:
                lengths[cur] = L
            cur = line[1:].split()[0]
            L = 0
        else:
            L += len(line)
    if cur is not None:
        lengths[cur] = L
    return lengths

te_len_map = read_fasta_lengths(te_fasta)

# load TE domain lookup
cdd_ids = set()
pfam_names = set()
with open(lookup_file) as f:
    _ = f.readline()  # skip header
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        _, db, name, domid = parts[:4]
        if db == "cdd":
            cdd_ids.add(domid)
        elif db == "pfam":
            pfam_names.add(name.lower())

# parse transeq frame suffix (_1 to _6)
_frame_re = re.compile(r"_(\d)$")

def parse_transeq_frame(q):
    m = _frame_re.search(q)
    if m:
        return q[:-2], int(m.group(1))
    return q, None

def pfam_aa_to_te_nt(aa_from, aa_to, frame6, te_len):
    if frame6 is None:
        return None, None, None
    if frame6 <= 3:
        strand = "+"
        offset = frame6 - 1
        nt_start = offset + (aa_from - 1) * 3 + 1
        nt_end = offset + aa_to * 3
    else:
        strand = "-"
        offset = frame6 - 4
        nt_end = te_len - offset - (aa_from - 1) * 3
        nt_start = te_len - offset - aa_to * 3 + 1
    if nt_start > nt_end:
        nt_start, nt_end = nt_end, nt_start
    nt_start = max(1, nt_start)
    nt_end = min(te_len, nt_end)
    return nt_start, nt_end, strand

hits = defaultdict(list)

# parse CDD hits (uses original TE IDs)
if cdd_file.exists():
    for line in cdd_file.open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 15:
            continue
        te_id = f[0]
        cdd_id = f[1].split("|")[-1]
        if cdd_id not in cdd_ids:
            continue
        qstart, qend = int(f[6]), int(f[7])
        sstart, send = int(f[8]), int(f[9])
        qframe = int(f[14])
        hits[te_id].append({
            "source": "CDD",
            "domain": cdd_id,
            "te_nt_start": min(qstart, qend),
            "te_nt_end": max(qstart, qend),
            "strand": "+" if qframe > 0 else "-",
            "frame": abs(qframe),
            "te_len": int(f[12]),
            "aa_start": ".",
            "aa_end": ".",
            "profile_start": min(sstart, send),
            "profile_end": max(sstart, send),
            "evalue": float(f[10]),
            "hit_q_full": te_id
        })

# parse Pfam hits (uses short IDs, need to map back)
if pfam_file.exists():
    for line in pfam_file.open():
        if line.startswith("#"):
            continue
        f = line.strip().split()
        if len(f) < 23:
            continue
        target = f[0]
        dom = f[3]
        if dom.lower() not in pfam_names:
            continue
        
        short_id, frame6 = parse_transeq_frame(target)
        full_te_id = short_to_full.get(short_id)
        if full_te_id is None:
            continue
        
        te_len = te_len_map.get(full_te_id)
        if te_len is None:
            continue
        
        hmm_from, hmm_to = int(f[15]), int(f[16])
        ali_from, ali_to = int(f[17]), int(f[18])
        nt0, nt1, strand = pfam_aa_to_te_nt(ali_from, ali_to, frame6, te_len)
        
        hits[full_te_id].append({
            "source": "Pfam",
            "domain": dom,
            "te_nt_start": nt0,
            "te_nt_end": nt1,
            "strand": strand,
            "frame": frame6,
            "te_len": te_len,
            "aa_start": ali_from,
            "aa_end": ali_to,
            "profile_start": hmm_from,
            "profile_end": hmm_to,
            "evalue": float(f[6]),
            "hit_q_full": target
        })

# write TSV
tsv = filtered_dir / f"{strain}_filtered_tes.tsv"
with open(tsv, "w") as out:
    out.write("te_id\tsource\tdomain\tte_nt_start\tte_nt_end\tstrand\tframe\tte_len\taa_start\taa_end\tprofile_start\tprofile_end\tevalue\thit_q_full\n")
    for te_id, doms in sorted(hits.items()):
        for d in doms:
            out.write(
                f"{te_id}\t{d['source']}\t{d['domain']}\t{d['te_nt_start']}\t{d['te_nt_end']}\t{d['strand']}\t{d['frame']}\t{d['te_len']}\t"
                f"{d['aa_start']}\t{d['aa_end']}\t{d['profile_start']}\t{d['profile_end']}\t{d['evalue']}\t{d['hit_q_full']}\n"
            )

# write filtered fasta (only TEs with domain hits)
outfa = filtered_dir / f"{strain}_filtered_tes.fasta"
keep = set(hits.keys())
cur = None
buf = []

with te_fasta.open() as fin, outfa.open("w") as fout:
    for line in fin:
        if line.startswith(">"):
            if cur is not None and cur in keep:
                fout.write(f">{cur}\n{''.join(buf)}")
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if cur is not None and cur in keep:
        fout.write(f">{cur}\n{''.join(buf)}")

print(f"    TEs with domain hits: {len(keep)}")
EOF

done

echo ""
echo "STEP 5 complete."

conda deactivate



### ===========================================================================
### STEP 6: Deduplicate and combine TEs (single job mode only, or SUMMARIZE=true)
### ===========================================================================

if [ "${RUN_MODE}" = "single" ] || [ "${SUMMARIZE}" = "true" ]; then

echo ""
echo "### STEP 6: Deduplicating and combining TEs ###"
echo ""

conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"

export RESULTS

python3 << 'EOF'
import os
from pathlib import Path
from hashlib import md5

results_dir = Path(os.environ.get('RESULTS', '.'))
filtered_dir = results_dir / "filtered_tes"
dedup_dir = results_dir / "deduplicated"
summary_dir = results_dir / "summary"

dedup_dir.mkdir(exist_ok=True)
summary_dir.mkdir(exist_ok=True)

def read_fasta(filepath):
    """Read fasta file, yield (header, sequence) tuples."""
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if header is not None:
            yield header, "".join(seq_parts)

def write_fasta(filepath, records):
    """Write list of (header, sequence) to fasta file."""
    with open(filepath, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

# process each strain
stats = []
all_tes = []

for fasta in sorted(filtered_dir.glob("*_filtered_tes.fasta")):
    strain = fasta.stem.replace("_filtered_tes", "")
    
    sequences = list(read_fasta(fasta))
    before_count = len(sequences)
    
    # deduplicate by exact sequence hash
    seen_hashes = {}
    dedup_records = []
    
    for header, seq in sequences:
        seq_hash = md5(seq.encode()).hexdigest()
        if seq_hash not in seen_hashes:
            seen_hashes[seq_hash] = header
            dedup_records.append((header, seq))
    
    after_count = len(dedup_records)
    removed = before_count - after_count
    
    # write deduplicated fasta
    outfile = dedup_dir / f"{strain}_dedup.fasta"
    write_fasta(outfile, dedup_records)
    
    # add to combined list with strain prefix
    for header, seq in dedup_records:
        new_id = f"{strain}|{header}"
        all_tes.append((new_id, strain, header, seq))
    
    stats.append({
        "strain": strain,
        "before": before_count,
        "after": after_count,
        "removed": removed
    })
    
    print(f"  {strain}: {before_count} -> {after_count} TEs (removed {removed} identical)")

# write combined fasta
combined_file = results_dir / "all_tes_combined.fasta"
with open(combined_file, "w") as f:
    for new_id, strain, orig_id, seq in all_tes:
        f.write(f">{new_id}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

# write info table
info_file = results_dir / "all_tes_info.tsv"
with open(info_file, "w") as f:
    f.write("te_id\tstrain\toriginal_id\n")
    for new_id, strain, orig_id, seq in all_tes:
        f.write(f"{new_id}\t{strain}\t{orig_id}\n")

# write dedup summary
summary_file = summary_dir / "dedup_summary.tsv"
with open(summary_file, "w") as f:
    f.write("strain\tbefore_dedup\tafter_dedup\tremoved\n")
    total_before = 0
    total_after = 0
    total_removed = 0
    for s in stats:
        f.write(f"{s['strain']}\t{s['before']}\t{s['after']}\t{s['removed']}\n")
        total_before += s['before']
        total_after += s['after']
        total_removed += s['removed']
    f.write(f"TOTAL\t{total_before}\t{total_after}\t{total_removed}\n")

print(f"\nTotal: {total_before:,} -> {total_after:,} TEs")
print(f"Removed {total_removed:,} identical sequences ({total_removed/total_before*100:.1f}%)")
EOF

echo ""
echo "STEP 6 complete."

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
echo "  Header maps:          ${RESULTS}/header_maps/"
echo "  Translated seqs:      ${RESULTS}/translated/"
echo "  CDD hits:             ${RESULTS}/cdd_hits/"
echo "  Pfam hits:            ${RESULTS}/pfam_hits/"
echo "  Filtered TEs:         ${RESULTS}/filtered_tes/"
echo "  Deduplicated TEs:     ${RESULTS}/deduplicated/"
echo "  Combined FASTA:       ${RESULTS}/all_tes_combined.fasta"
echo "  TE info table:        ${RESULTS}/all_tes_info.tsv"
echo "  Dedup summary:        ${RESULTS}/summary/dedup_summary.tsv"
echo ""
echo "Usage notes:"
echo "  - For array jobs, run dedup/combine separately after all tasks complete:"
echo "    SUMMARIZE=true sbatch 05_te_classification.sh"
echo "  - Or run as single job to process all genomes sequentially"
echo ""