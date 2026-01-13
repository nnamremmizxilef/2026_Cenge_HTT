#!/bin/bash
#SBATCH --job-name=te_classify
#SBATCH --qos hprio
#SBATCH --account node
#SBATCH --partition node
#SBATCH --mail-user=felix.zimmermann@wsl.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-22%11
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=5-00:00:00
#SBATCH --output /dev/null
#SBATCH --error /dev/null

set -euo pipefail

ROOT=/storage/zimmermf/HTT
TE_RESULTS=${ROOT}/results/04_results_te_annotation
RESULTS=${ROOT}/results/05_results_te_classification
LOGS=${ROOT}/logs/05_logs_te_classification

CDD_DB=${ROOT}/databases/cdd/Cdd
PFAM_DB=${ROOT}/databases/pfam/Pfam-A.TEonly.hmm
TE_DOMAINS_LOOKUP=${ROOT}/data/te_domains_lookup.tsv

mkdir -p "${RESULTS}"/{translated,cdd_hits,pfam_hits,filtered_tes,header_maps}
mkdir -p "${LOGS}"

# pick genome by array index (0-based)
mapfile -t TE_FASTAS < <(ls -1 "${TE_RESULTS}/te_sequences/"*_te.fasta | sort)
TE_COUNT=${#TE_FASTAS[@]}

if [ "${SLURM_ARRAY_TASK_ID}" -ge "${TE_COUNT}" ]; then
  echo "Task ID ${SLURM_ARRAY_TASK_ID} exceeds TE fasta count ${TE_COUNT}, exiting."
  exit 0
fi

TE_FASTA="${TE_FASTAS[$SLURM_ARRAY_TASK_ID]}"
STRAIN_ID=$(basename "${TE_FASTA}" _te.fasta)

# logging per genome
JOBTAG="${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${STRAIN_ID}"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "Strain: ${STRAIN_ID}"
echo "TE fasta: ${TE_FASTA}"

eval "$(conda shell.bash hook)"

### STEP 0: Create simplified fasta with short IDs for transeq
# transeq truncates long headers, so we create a mapping
SIMPLIFIED_FASTA="${RESULTS}/translated/${STRAIN_ID}_te_simplified.fasta"
HEADER_MAP="${RESULTS}/header_maps/${STRAIN_ID}_header_map.tsv"

if [ ! -f "${HEADER_MAP}" ]; then
  echo "Creating simplified headers for transeq..."
  awk '
  BEGIN { n=0 }
  /^>/ {
    n++
    full = substr($0, 2)
    # remove everything after first space
    split(full, arr, " ")
    full = arr[1]
    short = "TE" n
    print short "\t" full > "'"${HEADER_MAP}"'"
    print ">" short
    next
  }
  { print }
  ' "${TE_FASTA}" > "${SIMPLIFIED_FASTA}"
  echo "Header map: ${HEADER_MAP}"
fi

### STEP 1: transeq on simplified fasta
conda activate HTT_emboss
TRANSLATED="${RESULTS}/translated/${STRAIN_ID}_te_translated.faa"
if [ ! -f "${TRANSLATED}" ]; then
  echo "Running transeq..."
  transeq -sequence "${SIMPLIFIED_FASTA}" -outseq "${TRANSLATED}" -frame 6 -clean
fi
conda deactivate

### STEP 2: CDD rpstblastn (on original fasta - rpstblastn handles translation)
conda activate HTT_blast
CDD_OUT="${RESULTS}/cdd_hits/${STRAIN_ID}_cdd_hits.tsv"
if [ ! -f "${CDD_OUT}" ]; then
  echo "Running rpstblastn..."
  rpstblastn \
    -query "${TE_FASTA}" \
    -db "${CDD_DB}" \
    -out "${CDD_OUT}" \
    -evalue 0.001 \
    -max_target_seqs 1 \
    -num_threads "${SLURM_CPUS_PER_TASK}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qframe"
fi
conda deactivate

### STEP 3: Pfam hmmsearch on translated sequences
conda activate HTT_hmmer
PFAM_OUT="${RESULTS}/pfam_hits/${STRAIN_ID}_pfam_hits.tbl"
PFAM_DOM="${RESULTS}/pfam_hits/${STRAIN_ID}_pfam_domtbl.tbl"
if [ ! -f "${PFAM_DOM}" ]; then
  echo "Running hmmsearch..."
  hmmsearch \
    --tblout "${PFAM_OUT}" \
    --domtblout "${PFAM_DOM}" \
    -E 0.001 \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    "${PFAM_DB}" \
    "${TRANSLATED}"
fi
conda deactivate

### STEP 4: filter + enriched TSV + filtered fasta
conda activate HTT_python

python3 - "${TE_FASTA}" "${CDD_OUT}" "${PFAM_DOM}" "${RESULTS}" "${TE_DOMAINS_LOOKUP}" "${STRAIN_ID}" "${HEADER_MAP}" << 'PY'
import sys
from pathlib import Path
from collections import defaultdict
import re

te_fasta = Path(sys.argv[1])
cdd_file = Path(sys.argv[2])
pfam_file = Path(sys.argv[3])
results_dir = Path(sys.argv[4])
lookup_file = Path(sys.argv[5])
strain = sys.argv[6]
header_map_file = Path(sys.argv[7])

filtered_dir = results_dir / "filtered_tes"
filtered_dir.mkdir(exist_ok=True)

# Load header mapping (short -> full)
short_to_full = {}
with open(header_map_file) as f:
    for line in f:
        short, full = line.rstrip("\n").split("\t")
        short_to_full[short] = full

# Read original fasta to get TE lengths
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

# Load TE domain lookup
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

# Parse transeq frame suffix (_1 to _6)
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

# Parse CDD hits (uses original TE IDs)
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

# Parse Pfam hits (uses short IDs, need to map back)
if pfam_file.exists():
    for line in pfam_file.open():
        if line.startswith("#"):
            continue
        f = line.strip().split()
        if len(f) < 23:
            continue
        # Column 0 = target name (our translated seq), Column 3 = query name (Pfam domain)
        target = f[0]
        dom = f[3]
        if dom.lower() not in pfam_names:
            continue
        
        # Parse short ID and frame from target (e.g., "TE1_1" -> "TE1", frame 1)
        short_id, frame6 = parse_transeq_frame(target)
        
        # Map short ID back to full TE ID
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

# Write TSV
tsv = filtered_dir / f"{strain}_filtered_tes.tsv"
with open(tsv, "w") as out:
    out.write("te_id\tsource\tdomain\tte_nt_start\tte_nt_end\tstrand\tframe\tte_len\taa_start\taa_end\tprofile_start\tprofile_end\tevalue\thit_q_full\n")
    for te_id, doms in sorted(hits.items()):
        for d in doms:
            out.write(
                f"{te_id}\t{d['source']}\t{d['domain']}\t{d['te_nt_start']}\t{d['te_nt_end']}\t{d['strand']}\t{d['frame']}\t{d['te_len']}\t"
                f"{d['aa_start']}\t{d['aa_end']}\t{d['profile_start']}\t{d['profile_end']}\t{d['evalue']}\t{d['hit_q_full']}\n"
            )

# Write filtered fasta (only TEs with domain hits)
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
    # Write last sequence
    if cur is not None and cur in keep:
        fout.write(f">{cur}\n{''.join(buf)}")

print(f"TEs with domain hits: {len(keep)}")
print(f"Output TSV: {tsv}")
print(f"Output FASTA: {outfa}")
PY

conda deactivate
echo "Done strain ${STRAIN_ID}"