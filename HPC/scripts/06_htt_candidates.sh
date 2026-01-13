#!/bin/bash
#
#SBATCH --job-name=htt_candidates
#SBATCH --qos normal
#SBATCH --account node
#SBATCH --partition bigmem
#SBATCH --mail-user=felix.zimmermann@wsl.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output /storage/zimmermf/HTT/logs/06_logs_htt_candidates/htt_candidates_%j.out
#SBATCH --error  /storage/zimmermf/HTT/logs/06_logs_htt_candidates/htt_candidates_%j.err

set -euo pipefail

###############################################################################
### PATH DEFINITIONS
###############################################################################
ROOT=/storage/zimmermf/HTT
TREE_RESULTS=${ROOT}/results/01_results_phylo_tree
KS_RESULTS=${ROOT}/results/03_results_ks_divergence
TE_CLASS_RESULTS=${ROOT}/results/05_results_te_classification
RESULTS=${ROOT}/results/06_results_htt_candidates
LOGS=${ROOT}/logs/06_logs_htt_candidates

# Input files
TREE_FILE=${TREE_RESULTS}/phylogeny_rooted.treefile
KS_SUMMARY=${KS_RESULTS}/ks_calculations/node_ks_summary.tsv
NODES_INFO=${KS_RESULTS}/node_analyses/nodes_info.json
TE_FASTA=${TE_CLASS_RESULTS}/all_tes_combined.fasta
TE_INFO=${TE_CLASS_RESULTS}/all_tes_info.tsv
TE_DOMAINS_DIR=${TE_CLASS_RESULTS}/filtered_tes
CLADE_ASSIGN=${ROOT}/data/clade_assignment.tsv

# Intermediate files
BLAST_FASTA=${RESULTS}/blast_db/te_db.fasta
ID_MAP=${RESULTS}/blast_db/te_id_map.tsv
BLAST_OUT=${RESULTS}/blast_hits/all_vs_all_blast.tsv
BLAST_GZ=${RESULTS}/blast_hits/all_vs_all_blast.tsv.gz
FILTERED_BLAST=${RESULTS}/blast_hits/filtered_blast_both_directions.tsv
SORTED_BLAST=${RESULTS}/blast_hits/filtered_blast_sorted.tsv
CHAINED_FILE=${RESULTS}/chained_hits/chained_hits.tsv
PAIRS_JSON=${RESULTS}/ks_calculations/te_pairs_for_ks.json
CODON_DIR=${RESULTS}/ks_calculations/codon_alignments
TMP_DIR=${RESULTS}/ks_calculations/tmp
MANIFEST=${TMP_DIR}/pair_manifest.tsv
TE_KS_OUT=${RESULTS}/ks_calculations/te_ks_values.tsv
TE_KS_BACKUP=${RESULTS}/ks_calculations/te_ks_values_raw_backup.tsv

### Parameters (as per Romeijn et al. 2025)
# BLAST
BLAST_IDENTITY=75
BLAST_LENGTH=300
BLAST_BITSCORE=200
BLAST_EVALUE=10
BLAST_MAX_TARGETS=100000
BLAST_DBSIZE=10000000000

# Chaining
MAX_GAP=600
MAX_OVERLAP=600
MIN_COVERAGE=60

# Cores
NCORES="${SLURM_CPUS_PER_TASK:-64}"

# Ks calculation settings (optimized for large datasets)
KS_BATCH_N=1000
KS_CHUNK_FILES=200000

###############################################################################
### SETUP
###############################################################################
mkdir -p "${RESULTS}"/{blast_db,blast_hits,chained_hits,ks_calculations/codon_alignments,ks_calculations/tmp,htt_candidates,summary}
mkdir -p "${LOGS}"

JOBTAG="${SLURM_JOB_ID}_htt_candidates"
exec > >(tee "${LOGS}/${JOBTAG}.log") 2>&1

echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

export LANG=C
export LC_ALL=C
unset LANGUAGE
export PYTHONWARNINGS="ignore::FutureWarning"

# Prevent thread oversubscription (critical for R parallel)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

eval "$(conda shell.bash hook)"

START_TIME=$(date +%s)

echo ""
echo "Running HTT detection pipeline on HYPERION - $(date)"
echo "Results directory: ${RESULTS}"
echo "Using ${NCORES} cores"
echo ""
echo "Parameters:"
echo "  BLAST: identity >= ${BLAST_IDENTITY}%, length >= ${BLAST_LENGTH} bp, bitscore >= ${BLAST_BITSCORE}"
echo "  Chaining: max gap/overlap = ${MAX_GAP}/${MAX_OVERLAP} bp, min coverage = ${MIN_COVERAGE}%"
echo ""

###############################################################################
### SANITY CHECKS
###############################################################################
for f in "${TE_FASTA}" "${TE_INFO}" "${KS_SUMMARY}" "${NODES_INFO}" "${CLADE_ASSIGN}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

TE_COUNT=$(grep -c "^>" "${TE_FASTA}")
echo "Found ${TE_COUNT} TEs in database"

###############################################################################
### STEP 1: Create simplified FASTA and BLAST database
###############################################################################
echo ""
echo "### STEP 1: Creating BLAST database ###"
echo ""

if [ -f "${RESULTS}/blast_db/te_db.ndb" ]; then
    echo "BLAST database exists, skipping..."
else
    conda activate HTT_python

    python3 - "${TE_FASTA}" "${BLAST_FASTA}" "${ID_MAP}" << 'EOF'
import sys
from pathlib import Path

te_fasta = Path(sys.argv[1])
out_fasta = Path(sys.argv[2])
map_file = Path(sys.argv[3])

out_fasta.parent.mkdir(parents=True, exist_ok=True)

n = 0
with open(te_fasta) as f_in, open(out_fasta, "w") as f_out, open(map_file, "w") as f_map:
    f_map.write("short_id\tfull_id\n")
    for line in f_in:
        if line.startswith(">"):
            n += 1
            full_id = line[1:].split()[0]
            short_id = f"TE{n:07d}"
            f_map.write(f"{short_id}\t{full_id}\n")
            f_out.write(f">{short_id}\n")
        else:
            f_out.write(line)

print(f"Created simplified FASTA with {n} sequences")
EOF

    conda deactivate

    conda activate HTT_blast
    makeblastdb \
        -in "${BLAST_FASTA}" \
        -dbtype nucl \
        -out "${RESULTS}/blast_db/te_db" \
        -parse_seqids
    conda deactivate
    
    echo "BLAST database created"
fi

echo "STEP 1 complete."

###############################################################################
### STEP 2: All-vs-all BLASTn search
###############################################################################
echo ""
echo "### STEP 2: Running all-vs-all BLASTn ###"
echo ""

if [ -f "${BLAST_GZ}" ]; then
    echo "BLAST output exists (gzipped), skipping..."
elif [ -f "${BLAST_OUT}" ]; then
    echo "BLAST output exists, compressing..."
    gzip "${BLAST_OUT}"
else
    conda activate HTT_blast
    
    blastn \
        -query "${BLAST_FASTA}" \
        -db "${RESULTS}/blast_db/te_db" \
        -out "${BLAST_OUT}" \
        -dbsize ${BLAST_DBSIZE} \
        -evalue ${BLAST_EVALUE} \
        -max_target_seqs ${BLAST_MAX_TARGETS} \
        -num_threads "${NCORES}" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

    conda deactivate
    
    echo "BLAST complete, compressing..."
    gzip "${BLAST_OUT}"
fi

echo "STEP 2 complete."

###############################################################################
### STEP 3: Filter BLAST hits (keeping both directions for chaining)
###############################################################################
echo ""
echo "### STEP 3: Filtering BLAST hits ###"
echo ""

if [ -f "${FILTERED_BLAST}" ] && [ -s "${FILTERED_BLAST}" ]; then
    echo "Filtered BLAST file exists, skipping..."
    echo "Lines: $(wc -l < ${FILTERED_BLAST})"
else
    conda activate HTT_python

    python3 - "${BLAST_GZ}" "${ID_MAP}" "${TE_INFO}" "${CLADE_ASSIGN}" "${FILTERED_BLAST}" \
             "${BLAST_IDENTITY}" "${BLAST_LENGTH}" "${BLAST_BITSCORE}" << 'PYEOF'
import sys
import gzip
from pathlib import Path

blast_file = Path(sys.argv[1])
id_map_file = Path(sys.argv[2])
te_info_file = Path(sys.argv[3])
clade_file = Path(sys.argv[4])
filtered_file = Path(sys.argv[5])
min_identity = float(sys.argv[6])
min_length = int(sys.argv[7])
min_bitscore = float(sys.argv[8])

# Load short -> full ID mapping
short_to_full = {}
with open(id_map_file) as f:
    _ = f.readline()
    for line in f:
        short, full = line.rstrip("\n").split("\t")
        short_to_full[short] = full
print(f"Loaded ID mapping for {len(short_to_full)} TEs")

# TE (full ID) -> genome
te_genome = {}
with open(te_info_file) as f:
    _ = f.readline()
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        te_genome[fields[0]] = fields[1]
print(f"Loaded genome info for {len(te_genome)} TEs")

# genome -> clade
genome_to_clade = {}
with open(clade_file) as f:
    header = f.readline().strip().split("\t")
    col = {h:i for i,h in enumerate(header)}
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) < 2:
            continue
        sample = fields[col.get("sample", 0)]
        clade = fields[col.get("clade", 1)]
        genome_to_clade[sample] = clade
print(f"Loaded clade assignments for {len(genome_to_clade)} genomes")

def get_clade(genome):
    return genome_to_clade.get(genome, genome)

hits_total = 0
hits_passed = 0

with gzip.open(blast_file, 'rt') as f_in, open(filtered_file, "w") as f_out:
    f_out.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend"
                "\tsstart\tsend\tevalue\tbitscore\tqlen\tslen"
                "\tqgenome\tsgenome\tqclade\tsclade\n")

    for line in f_in:
        hits_total += 1
        
        if hits_total % 100000000 == 0:
            print(f"  Processed {hits_total:,} hits, passed {hits_passed:,}", flush=True)
        
        fields = line.strip().split("\t")
        if len(fields) < 14:
            continue

        qshort, sshort = fields[0], fields[1]
        
        # Skip self-hits
        if qshort == sshort:
            continue

        pident = float(fields[2])
        length = int(fields[3])
        bitscore = float(fields[11])

        # Apply filters
        if pident < min_identity or length < min_length or bitscore < min_bitscore:
            continue

        qseqid = short_to_full.get(qshort, qshort)
        sseqid = short_to_full.get(sshort, sshort)

        qgenome = te_genome.get(qseqid, "unknown")
        sgenome = te_genome.get(sseqid, "unknown")
        qclade = get_clade(qgenome)
        sclade = get_clade(sgenome)

        # Only keep hits between different clades
        if qclade == sclade:
            continue

        hits_passed += 1
        f_out.write(f"{qseqid}\t{sseqid}\t{fields[2]}\t{fields[3]}\t{fields[4]}\t{fields[5]}\t"
                    f"{fields[6]}\t{fields[7]}\t{fields[8]}\t{fields[9]}\t{fields[10]}\t{fields[11]}\t"
                    f"{fields[12]}\t{fields[13]}\t{qgenome}\t{sgenome}\t{qclade}\t{sclade}\n")

print(f"")
print(f"Total BLAST hits: {hits_total:,}")
print(f"Filtered hits (both directions): {hits_passed:,}")
PYEOF

    conda deactivate
fi

echo "STEP 3 complete."

###############################################################################
### STEP 4: Sort and chain BLAST hits, then apply one-direction filter
###############################################################################
echo ""
echo "### STEP 4: Sorting and chaining BLAST hits ###"
echo ""

# Step 4a: Sort by qseqid, sseqid
if [ -f "${SORTED_BLAST}" ] && [ -s "${SORTED_BLAST}" ]; then
    echo "Sorted BLAST file exists, skipping sort..."
else
    echo "Sorting filtered BLAST hits..."
    
    # Keep header, sort rest
    head -1 "${FILTERED_BLAST}" > "${SORTED_BLAST}"
    tail -n +2 "${FILTERED_BLAST}" | sort -t$'\t' -k1,1 -k2,2 -S 50G --parallel=${NCORES} >> "${SORTED_BLAST}"
    
    echo "Sorting complete"
fi

# Step 4b: Chain and apply one-direction filter
if [ -f "${CHAINED_FILE}" ] && [ -s "${CHAINED_FILE}" ]; then
    echo "Chained hits file exists, skipping..."
    echo "Lines: $(wc -l < ${CHAINED_FILE})"
else
    echo "Chaining hits..."
    
    conda activate HTT_python

    python3 - "${SORTED_BLAST}" "${CHAINED_FILE}" "${MAX_GAP}" "${MAX_OVERLAP}" "${MIN_COVERAGE}" << 'PYEOF'
import sys
from pathlib import Path

sorted_file = Path(sys.argv[1])
chained_file = Path(sys.argv[2])
max_gap = int(sys.argv[3])
max_overlap = int(sys.argv[4])
min_coverage = float(sys.argv[5])

def chain_hsps(hsps, max_gap, max_overlap):
    """Chain HSPs sorted by query start position."""
    if not hsps:
        return []
    
    hsps_sorted = sorted(hsps, key=lambda x: x["qstart"])
    chains = []
    
    for hsp in hsps_sorted:
        added = False
        for chain in chains:
            last = chain[-1]
            q_gap = hsp["qstart"] - last["qend"]
            s_gap = hsp["sstart"] - last["send"]
            
            if -max_overlap <= q_gap <= max_gap and -max_overlap <= s_gap <= max_gap:
                chain.append(hsp)
                added = True
                break
        
        if not added:
            chains.append([hsp])
    
    return chains

def merge_chain(chain):
    """Merge a chain of HSPs into a single hit."""
    qstart = min(h["qstart"] for h in chain)
    qend = max(h["qend"] for h in chain)
    sstart = min(h["sstart"] for h in chain)
    send = max(h["send"] for h in chain)
    
    return {
        "qstart": qstart, "qend": qend,
        "sstart": sstart, "send": send,
        "qlen": chain[0]["qlen"], "slen": chain[0]["slen"],
        "bitscore": sum(h["bitscore"] for h in chain),
        "pident": sum(h["pident"] * h["length"] for h in chain) / sum(h["length"] for h in chain),
        "qgenome": chain[0]["qgenome"], "sgenome": chain[0]["sgenome"]
    }

def chain_and_filter(qseqid, sseqid, hsps):
    """Chain HSPs and filter by coverage."""
    chains = chain_hsps(hsps, max_gap, max_overlap)
    
    results = []
    for chain in chains:
        merged = merge_chain(chain)
        
        q_cov = 100.0 * (merged["qend"] - merged["qstart"] + 1) / merged["qlen"]
        s_cov = 100.0 * (merged["send"] - merged["sstart"] + 1) / merged["slen"]
        
        if q_cov >= min_coverage and s_cov >= min_coverage:
            merged["q_coverage"] = q_cov
            merged["s_coverage"] = s_cov
            results.append(merged)
    
    return results

# Process sorted file
total_lines = 0
total_pairs = 0
total_chains_before_dedup = 0
total_chains_after_dedup = 0
skipped_lines = 0
seen_pairs = set()

current_pair = None
current_hits = []

with open(sorted_file) as f_in, open(chained_file, "w") as f_out:
    header = f_in.readline()
    f_out.write("qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tqlen\tslen\t"
                "q_coverage\ts_coverage\tbitscore\tpident\tqgenome\tsgenome\n")
    
    for line in f_in:
        total_lines += 1
        fields = line.strip().split("\t")
        if len(fields) < 16:
            skipped_lines += 1
            continue
        
        pair = (fields[0], fields[1])
        
        if pair != current_pair:
            # Process previous pair
            if current_hits:
                results = chain_and_filter(current_pair[0], current_pair[1], current_hits)
                total_chains_before_dedup += len(results)
                
                for r in results:
                    q, s = current_pair[0], current_pair[1]
                    canonical_pair = (min(q, s), max(q, s))
                    if canonical_pair in seen_pairs:
                        continue
                    seen_pairs.add(canonical_pair)
                    
                    # Write in canonical order (alphabetically smaller first)
                    if q > s:
                        f_out.write(
                            f"{s}\t{q}\t{r['sstart']}\t{r['send']}\t"
                            f"{r['qstart']}\t{r['qend']}\t{r['slen']}\t{r['qlen']}\t"
                            f"{r['s_coverage']:.3f}\t{r['q_coverage']:.3f}\t"
                            f"{r['bitscore']:.1f}\t{r['pident']:.2f}\t"
                            f"{r['sgenome']}\t{r['qgenome']}\n"
                        )
                    else:
                        f_out.write(
                            f"{q}\t{s}\t{r['qstart']}\t{r['qend']}\t"
                            f"{r['sstart']}\t{r['send']}\t{r['qlen']}\t{r['slen']}\t"
                            f"{r['q_coverage']:.3f}\t{r['s_coverage']:.3f}\t"
                            f"{r['bitscore']:.1f}\t{r['pident']:.2f}\t"
                            f"{r['qgenome']}\t{r['sgenome']}\n"
                        )
                    total_chains_after_dedup += 1
                
                total_pairs += 1
            
            current_pair = pair
            current_hits = []
        
        try:
            current_hits.append({
                "qstart": int(fields[6]), "qend": int(fields[7]),
                "sstart": int(fields[8]), "send": int(fields[9]),
                "qlen": int(fields[12]), "slen": int(fields[13]),
                "bitscore": float(fields[11]), "pident": float(fields[2]),
                "length": int(fields[3]),
                "qgenome": fields[14], "sgenome": fields[15]
            })
        except (ValueError, IndexError):
            skipped_lines += 1
            continue
        
        if total_lines % 50000000 == 0:
            print(f"  Processed {total_lines:,} lines, {total_pairs:,} pairs, "
                  f"{total_chains_after_dedup:,} chains", flush=True)
    
    # Process last pair
    if current_hits:
        results = chain_and_filter(current_pair[0], current_pair[1], current_hits)
        total_chains_before_dedup += len(results)
        
        for r in results:
            q, s = current_pair[0], current_pair[1]
            canonical_pair = (min(q, s), max(q, s))
            if canonical_pair in seen_pairs:
                continue
            seen_pairs.add(canonical_pair)
            
            if q > s:
                f_out.write(
                    f"{s}\t{q}\t{r['sstart']}\t{r['send']}\t"
                    f"{r['qstart']}\t{r['qend']}\t{r['slen']}\t{r['qlen']}\t"
                    f"{r['s_coverage']:.3f}\t{r['q_coverage']:.3f}\t"
                    f"{r['bitscore']:.1f}\t{r['pident']:.2f}\t"
                    f"{r['sgenome']}\t{r['qgenome']}\n"
                )
            else:
                f_out.write(
                    f"{q}\t{s}\t{r['qstart']}\t{r['qend']}\t"
                    f"{r['sstart']}\t{r['send']}\t{r['qlen']}\t{r['slen']}\t"
                    f"{r['q_coverage']:.3f}\t{r['s_coverage']:.3f}\t"
                    f"{r['bitscore']:.1f}\t{r['pident']:.2f}\t"
                    f"{r['qgenome']}\t{r['sgenome']}\n"
                )
            total_chains_after_dedup += 1
        
        total_pairs += 1

print(f"")
print(f"Total lines processed: {total_lines:,}")
print(f"Total TE pairs (both directions): {total_pairs:,}")
print(f"Chained hits before one-direction filter: {total_chains_before_dedup:,}")
print(f"Chained hits after one-direction filter: {total_chains_after_dedup:,}")
PYEOF

    conda deactivate
fi

echo "STEP 4 complete."

###############################################################################
### STEP 5: Select TE pairs with shared domains for Ks calculation
###############################################################################
echo ""
echo "### STEP 5: Extracting TE pairs with shared domains ###"
echo ""

if [ -f "${PAIRS_JSON}" ] && [ -s "${PAIRS_JSON}" ]; then
    echo "Pairs JSON exists, skipping..."
else
    conda activate HTT_python

    python3 - "${CHAINED_FILE}" "${TE_DOMAINS_DIR}" "${PAIRS_JSON}" << 'PYEOF'
import sys
from pathlib import Path
from collections import defaultdict
import json

chained_file = Path(sys.argv[1])
te_domains_dir = Path(sys.argv[2])
pairs_file = Path(sys.argv[3])

# Load domain annotations
te_domains = defaultdict(list)

for domain_file in te_domains_dir.glob("*_filtered_tes.tsv"):
    genome_name = domain_file.stem.replace("_filtered_tes", "")
    with open(domain_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h:i for i,h in enumerate(header)}
        required = ["te_id","source","domain","te_nt_start","te_nt_end","profile_start","profile_end","strand"]
        for r in required:
            if r not in idx:
                raise SystemExit(f"ERROR: {domain_file} missing required column: {r}")

        for line in f:
            fields = line.rstrip("\n").split("\t")
            te_id = fields[idx["te_id"]]
            full_te_id = f"{genome_name}|{te_id}"

            def _int(name):
                v = fields[idx[name]]
                return None if v in ("", "None", ".", "NA") else int(v)

            te_nt_start = _int("te_nt_start")
            te_nt_end   = _int("te_nt_end")
            prof_start  = _int("profile_start")
            prof_end    = _int("profile_end")
            strand      = fields[idx["strand"]]
            source      = fields[idx["source"]]
            dom         = fields[idx["domain"]]

            if te_nt_start is None or te_nt_end is None or prof_start is None or prof_end is None:
                continue

            te_domains[full_te_id].append({
                "source": source,
                "domain": dom,
                "te_start": min(te_nt_start, te_nt_end),
                "te_end": max(te_nt_start, te_nt_end),
                "strand": strand,
                "profile_start": prof_start,
                "profile_end": prof_end
            })

print(f"Loaded domain annotations for {len(te_domains)} TEs")

# Load chained hits
chained_hits = []
with open(chained_file) as f:
    header = f.readline()
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) < 14:
            continue
        chained_hits.append({
            "qseqid": fields[0], "sseqid": fields[1],
            "qstart": int(fields[2]), "qend": int(fields[3]),
            "sstart": int(fields[4]), "send": int(fields[5]),
            "qgenome": fields[12], "sgenome": fields[13]
        })
print(f"Loaded {len(chained_hits)} chained hits")

def inside(a_start, a_end, b_start, b_end):
    return a_start >= b_start and a_end <= b_end

# Find pairs with shared domains
pairs_for_ks = []

for hit in chained_hits:
    qseqid = hit["qseqid"]
    sseqid = hit["sseqid"]

    q_dom = te_domains.get(qseqid, [])
    s_dom = te_domains.get(sseqid, [])
    if not q_dom or not s_dom:
        continue

    qh0, qh1 = min(hit["qstart"], hit["qend"]), max(hit["qstart"], hit["qend"])
    sh0, sh1 = min(hit["sstart"], hit["send"]), max(hit["sstart"], hit["send"])

    shared = []
    for qd in q_dom:
        if not inside(qd["te_start"], qd["te_end"], qh0, qh1):
            continue
        q_rel_start = qd["te_start"] - qh0
        q_rel_end   = qd["te_end"]   - qh0

        for sd in s_dom:
            # Must match: source, domain, profile coordinates
            if qd["source"] != sd["source"]:
                continue
            if qd["domain"] != sd["domain"]:
                continue
            if qd["profile_start"] != sd["profile_start"] or qd["profile_end"] != sd["profile_end"]:
                continue
            if not inside(sd["te_start"], sd["te_end"], sh0, sh1):
                continue

            s_rel_start = sd["te_start"] - sh0
            s_rel_end   = sd["te_end"]   - sh0

            # Must be at same relative position in alignment
            if q_rel_start != s_rel_start or q_rel_end != s_rel_end:
                continue

            shared.append({
                "source": qd["source"],
                "domain": qd["domain"],
                "profile_pos": (qd["profile_start"], qd["profile_end"]),
                "q_nt": (qd["te_start"], qd["te_end"], qd["strand"]),
                "s_nt": (sd["te_start"], sd["te_end"], sd["strand"]),
                "rel": (q_rel_start, q_rel_end)
            })

    if not shared:
        continue

    # Remove overlapping domains
    shared.sort(key=lambda x: x["rel"][0])
    non_overlapping = []
    last_end = -1
    total_bp = 0
    for d in shared:
        start, end = d["rel"]
        if start <= last_end:
            continue
        non_overlapping.append(d)
        total_bp += (d["q_nt"][1] - d["q_nt"][0] + 1)
        last_end = end

    # Require >= 300 bp of shared domains
    if total_bp >= 300:
        pairs_for_ks.append({
            "qseqid": qseqid,
            "sseqid": sseqid,
            "qgenome": hit["qgenome"],
            "sgenome": hit["sgenome"],
            "domains": non_overlapping,
            "total_domain_bp": total_bp
        })

print(f"TE pairs with shared domains for Ks calculation: {len(pairs_for_ks)}")

with open(pairs_file, "w") as out:
    json.dump(pairs_for_ks, out, indent=2)
PYEOF

    conda deactivate
fi

echo "STEP 5 complete."

###############################################################################
### STEP 6: Calculate Ks for TE pairs
###############################################################################
echo ""
echo "### STEP 6: Calculating TE Ks ###"
echo ""

# Step 6a: Prepare FASTA files for alignment
if [ -f "${MANIFEST}" ] && [ -s "${MANIFEST}" ]; then
    echo "Manifest exists, skipping preparation..."
    PAIR_COUNT=$(tail -n +2 "${MANIFEST}" | wc -l)
    echo "Found ${PAIR_COUNT} pairs in manifest"
else
    echo "Preparing FASTA files for alignment..."
    
    conda activate HTT_python

    python3 - "${PAIRS_JSON}" "${TE_FASTA}" "${TMP_DIR}" "${CODON_DIR}" << 'PYEOF'
import json, sys
from pathlib import Path

pairs_file = Path(sys.argv[1])
te_fasta = Path(sys.argv[2])
tmp_dir = Path(sys.argv[3])
codon_dir = Path(sys.argv[4])

pairs = json.loads(pairs_file.read_text())
if not pairs:
    print("No pairs for Ks; exiting.")
    sys.exit(0)

print(f"Loaded {len(pairs)} pairs from JSON")

# Load sequences
seqs = {}
cur = None
buf = []
with open(te_fasta) as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf).upper()
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf).upper()

print(f"Loaded {len(seqs)} TE sequences")

_comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s):
    return s.translate(_comp)[::-1]

CODON_TABLE = {
"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
"TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
"CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
"CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
"AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
"GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
"GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

def translate_nt(nt):
    return "".join(CODON_TABLE.get(nt[i:i+3], "X") for i in range(0, len(nt) - 2, 3))

def extract_domain(seq, start, end, strand):
    frag = seq[start-1:end]
    if strand == "-":
        frag = revcomp(frag)
    L = (len(frag) // 3) * 3
    return frag[:L]

ok = 0
skipped = {"seq_not_found": 0, "bad_frag": 0, "too_short": 0}

manifest = tmp_dir / "pair_manifest.tsv"
with open(manifest, "w") as man:
    man.write("pair_id\tqseqid\tsseqid\tnuc_fa\taa_fa\taa_aln\tcodon_aln\n")

    for i, pair in enumerate(pairs, start=1):
        qid, sid = pair["qseqid"], pair["sseqid"]
        doms = pair["domains"]

        if qid not in seqs or sid not in seqs:
            skipped["seq_not_found"] += 1
            continue

        qseq, sseq = seqs[qid], seqs[sid]
        q_parts, s_parts = [], []
        bad = False

        for d in doms:
            qs, qe, qstr = d["q_nt"]
            ss, se, sstr = d["s_nt"]
            qfrag = extract_domain(qseq, qs, qe, qstr)
            sfrag = extract_domain(sseq, ss, se, sstr)
            if len(qfrag) < 3 or len(sfrag) < 3:
                bad = True
                break
            m = (min(len(qfrag), len(sfrag)) // 3) * 3
            if m < 3:
                bad = True
                break
            q_parts.append(qfrag[:m])
            s_parts.append(sfrag[:m])

        if bad:
            skipped["bad_frag"] += 1
            continue

        q_concat, s_concat = "".join(q_parts), "".join(s_parts)
        if len(q_concat) < 300 or len(s_concat) < 300:
            skipped["too_short"] += 1
            continue

        q_aa = translate_nt(q_concat)
        s_aa = translate_nt(s_concat)
        
        # Replace stop codons with X (paper does not filter stop codons)
        q_aa = q_aa.replace("*", "X")
        s_aa = s_aa.replace("*", "X")

        pair_id = f"pair_{i:09d}"
        nuc_fa = tmp_dir / f"{pair_id}.nuc.fasta"
        aa_fa = tmp_dir / f"{pair_id}.aa.fasta"
        aa_aln = tmp_dir / f"{pair_id}.aa.aln.fasta"
        cod_aln = codon_dir / f"{pair_id}.codon_aln.fasta"

        with open(nuc_fa, "w") as out:
            out.write(f">query\n{q_concat}\n>subject\n{s_concat}\n")
        with open(aa_fa, "w") as out:
            out.write(f">query\n{q_aa}\n>subject\n{s_aa}\n")

        man.write(f"{pair_id}\t{qid}\t{sid}\t{nuc_fa}\t{aa_fa}\t{aa_aln}\t{cod_aln}\n")
        ok += 1

        if ok % 500000 == 0:
            print(f"  Prepared {ok} pairs...", flush=True)

print(f"Wrote {ok} pair inputs to manifest")
print(f"Skipped: {skipped}")
PYEOF

    conda deactivate
    
    PAIR_COUNT=$(tail -n +2 "${MANIFEST}" | wc -l)
    echo "Total pairs to process: ${PAIR_COUNT}"
fi

if [ "${PAIR_COUNT:-0}" -eq 0 ]; then
    echo "No pairs to process"
    exit 0
fi

# Step 6b: MAFFT alignment (parallel)
echo ""
echo "Running MAFFT (${NCORES} parallel jobs) at $(date)..."
conda activate HTT_mafft

tail -n +2 "${MANIFEST}" | \
    parallel -j ${NCORES} --colsep '\t' \
    '[ -f {7} ] || mafft --auto --quiet {5} > {6} 2>/dev/null'

echo "MAFFT complete at $(date)"
conda deactivate

# Step 6c: PAL2NAL (parallel)
echo ""
echo "Running PAL2NAL (${NCORES} parallel jobs) at $(date)..."
conda activate HTT_pal2nal

tail -n +2 "${MANIFEST}" | \
    parallel -j ${NCORES} --colsep '\t' \
    '[ -f {7} ] || [ ! -s {6} ] || pal2nal.pl {6} {4} -output fasta -nogap > {7} 2>/dev/null'

# Remove empty codon alignments
find "${CODON_DIR}" -name "*.fasta" -empty -delete 2>/dev/null || true

echo "PAL2NAL complete at $(date)"
conda deactivate

# Step 6d: Ks calculation using local SSD and GNU parallel (memory-efficient)
echo ""
echo "Calculating Ks at $(date)..."
conda activate HTT_r

# Setup local SSD for fast I/O
LOCAL_BASE="${TMPDIR:-/tmp}/htt_ks_${SLURM_JOB_ID}"
LOCAL_IN="${LOCAL_BASE}/in"
LOCAL_OUT="${LOCAL_BASE}/out"
mkdir -p "${LOCAL_IN}" "${LOCAL_OUT}"

# R batch script for efficient Ks calculation
cat > /tmp/calc_ks_batch.R << 'RSCRIPT'
library(seqinr)
files <- commandArgs(trailingOnly=TRUE)
for (f in files) {
  tryCatch({
    aln <- read.alignment(f, format="fasta")
    if (length(aln$seq) == 2) {
      ksres <- kaks(aln)
      pair_id <- sub("\\.codon_aln\\.fasta$", "", basename(f))
      cat(pair_id, ksres$ks, ksres$ka, "\n", sep="\t")
    }
  }, error=function(e) {})
}
RSCRIPT

# Get file list from manifest (avoids "argument list too long" error)
echo "Extracting file paths from manifest..."
FILELIST="${LOCAL_BASE}/codon_files.list"
cut -f7 "${MANIFEST}" | tail -n +2 > "${FILELIST}"
TOTAL=$(wc -l < "${FILELIST}")
echo "Total files: ${TOTAL}"

# Initialize output
echo -e "pair_id\tks\tka" > "${TE_KS_OUT}"

# Process in chunks to manage memory
START=1
CHUNK_NUM=0
while [ "${START}" -le "${TOTAL}" ]; do
    CHUNK_NUM=$((CHUNK_NUM + 1))
    END=$((START + KS_CHUNK_FILES - 1))
    if [ "${END}" -gt "${TOTAL}" ]; then END="${TOTAL}"; fi

    echo ""
    echo "=== Chunk ${CHUNK_NUM}: files ${START}-${END} / ${TOTAL} at $(date) ==="

    # Clean local dirs using find (avoids "argument list too long")
    find "${LOCAL_IN}" -type f -delete 2>/dev/null || true
    find "${LOCAL_OUT}" -type f -delete 2>/dev/null || true

    # Extract this chunk's file paths
    CHUNK_LIST="${LOCAL_BASE}/chunk.list"
    sed -n "${START},${END}p" "${FILELIST}" > "${CHUNK_LIST}"
    CHUNK_SIZE=$(wc -l < "${CHUNK_LIST}")
    echo "Chunk has ${CHUNK_SIZE} files"

    # Stage to local SSD
    echo "Staging to local SSD..."
    xargs -a "${CHUNK_LIST}" -P ${NCORES} -I{} cp -n "{}" "${LOCAL_IN}/" 2>/dev/null || true

    STAGED=$(find "${LOCAL_IN}" -type f | wc -l)
    echo "Staged ${STAGED} files to local SSD"

    # Run Ks calculation using find (not ls, avoids argument list errors)
    echo "Calculating Ks..."
    CHUNK_TSV="${LOCAL_OUT}/chunk_ks.tsv"
    
    find "${LOCAL_IN}" -type f -name "*.fasta" | \
        parallel -j ${NCORES} --pipe -N ${KS_BATCH_N} \
        'xargs Rscript /tmp/calc_ks_batch.R' > "${CHUNK_TSV}" 2>/dev/null || true

    # Append to main output
    CHUNK_RESULTS=$(wc -l < "${CHUNK_TSV}" || echo 0)
    echo "Chunk produced ${CHUNK_RESULTS} Ks values"
    cat "${CHUNK_TSV}" >> "${TE_KS_OUT}"

    TOTAL_SO_FAR=$(($(wc -l < "${TE_KS_OUT}") - 1))
    echo "Total Ks values so far: ${TOTAL_SO_FAR}"

    START=$((END + 1))
done

echo ""
echo "Ks calculation complete at $(date)"
FINAL_COUNT=$(($(wc -l < "${TE_KS_OUT}") - 1))
echo "Final Ks values: ${FINAL_COUNT}"

# Cleanup local SSD
rm -rf "${LOCAL_BASE}"

conda deactivate

# Create backup before modifying
echo ""
echo "Creating backup of raw Ks output..."
cp "${TE_KS_OUT}" "${TE_KS_BACKUP}"
echo "Backup saved to: ${TE_KS_BACKUP}"

# Clean trailing whitespace (fix for potential parsing issues)
sed -i 's/[[:space:]]*$//' "${TE_KS_OUT}"

# Add qseqid/sseqid from manifest
echo ""
echo "Adding sequence IDs to Ks output..."
conda activate HTT_python

python3 - << 'PYEOF'
import pandas as pd
from pathlib import Path

ks_file = Path("/storage/zimmermf/HTT/results/06_results_htt_candidates/ks_calculations/te_ks_values.tsv")
manifest = Path("/storage/zimmermf/HTT/results/06_results_htt_candidates/ks_calculations/tmp/pair_manifest.tsv")

# Read with explicit whitespace handling
ks = pd.read_csv(ks_file, sep="\t", dtype=str, skipinitialspace=True)
ks.columns = ks.columns.str.strip()
for col in ks.columns:
    if ks[col].dtype == object:
        ks[col] = ks[col].str.strip()

ks["ks"] = pd.to_numeric(ks["ks"], errors="coerce")
ks["ka"] = pd.to_numeric(ks["ka"], errors="coerce")

print(f"Loaded {len(ks)} Ks values")

# Read manifest
man = pd.read_csv(manifest, sep="\t", usecols=["pair_id", "qseqid", "sseqid"], dtype=str)
man["pair_id"] = man["pair_id"].str.strip()
man["qseqid"] = man["qseqid"].str.strip()
man["sseqid"] = man["sseqid"].str.strip()

print(f"Loaded {len(man)} manifest rows")

# Check match
common = set(ks["pair_id"]) & set(man["pair_id"])
print(f"Matching pair_ids: {len(common)}")

if len(common) == 0:
    raise ValueError("No matching pair_ids between Ks output and manifest!")

# Merge
merged = ks.merge(man, on="pair_id", how="left")
merged = merged[["pair_id", "qseqid", "sseqid", "ks", "ka"]]

matched = merged["qseqid"].notna().sum()
print(f"Rows with qseqid after merge: {matched} / {len(merged)}")

if matched == 0:
    raise ValueError("Merge failed - no rows have qseqid!")

merged.to_csv(ks_file, sep="\t", index=False)
print(f"Wrote {len(merged)} rows to {ks_file}")
PYEOF

conda deactivate

echo ""
echo "STEP 6 complete at $(date)"

###############################################################################
### STEP 7: Identify HTT candidates
###############################################################################
echo ""
echo "### STEP 7: Identifying HTT candidates ###"
echo ""

conda activate HTT_python

python3 - << 'PYEOF'
import json
from pathlib import Path
import pandas as pd
import numpy as np

RESULTS = Path("/storage/zimmermf/HTT/results/06_results_htt_candidates")

te_ks_file = RESULTS / "ks_calculations" / "te_ks_values.tsv"
pairs_json = RESULTS / "ks_calculations" / "te_pairs_for_ks.json"
busco_ks_file = Path("/storage/zimmermf/HTT/results/03_results_ks_divergence/ks_calculations/node_ks_summary.tsv")
nodes_file = Path("/storage/zimmermf/HTT/results/03_results_ks_divergence/node_analyses/nodes_info.json")
clade_file = Path("/storage/zimmermf/HTT/data/clade_assignment.tsv")

htt_dir = RESULTS / "htt_candidates"
summary_dir = RESULTS / "summary"
htt_dir.mkdir(exist_ok=True)
summary_dir.mkdir(exist_ok=True)

def get_short_name(full_name):
    """Extract 'Cenge1005' from 'Cenge1005_1_AssemblyScaffolds_2024-05-13'"""
    return full_name.split("_")[0]

# Load Ks values
print("Loading Ks values...")
te_ks = pd.read_csv(te_ks_file, sep="\t")
print(f"Loaded {len(te_ks)} TE Ks values")

# Build pair_meta from JSON (memory-efficient)
print("Loading pair metadata...")
pair_meta = {}
with open(pairs_json) as f:
    pairs = json.load(f)
for p in pairs:
    key = p["qseqid"] + "||" + p["sseqid"]
    pair_meta[key] = {"qgenome": p["qgenome"], "sgenome": p["sgenome"]}
del pairs  # Free memory
print(f"Loaded {len(pair_meta)} pair metadata entries")

# Load BUSCO thresholds
busco = pd.read_csv(busco_ks_file, sep="\t")
if "pass_0.01" in busco.columns:
    busco = busco[busco["pass_0.01"] == True]
busco_map = dict(zip(busco["node"], busco["ks_0.5_quantile"]))
print(f"Loaded {len(busco_map)} BUSCO Ks thresholds")
print(f"Threshold range: {min(busco_map.values()):.4f} - {max(busco_map.values()):.4f}")

# Load clade assignments
clade_map = {}
with open(clade_file) as f:
    header = f.readline().strip().split("\t")
    col = {h:i for i,h in enumerate(header)}
    for line in f:
        fields = line.strip().split("\t")
        clade_map[fields[col["sample"]]] = fields[col["clade"]]
print(f"Loaded {len(clade_map)} clade assignments")

# Build node -> clades mapping
nodes = json.loads(nodes_file.read_text())
node_to_clades = {}
node_size = {}
for n in nodes:
    clades = set()
    for leaf in n["all_leaves"]:
        short_name = get_short_name(leaf)
        clades.add(clade_map.get(short_name, short_name))
    node_to_clades[n["node_name"]] = clades
    node_size[n["node_name"]] = len(clades)
print(f"Built clade sets for {len(node_to_clades)} nodes")

def find_mrca(c1, c2):
    hits = [n for n, cs in node_to_clades.items() if c1 in cs and c2 in cs]
    return min(hits, key=lambda n: node_size[n]) if hits else None

# Classify HTT candidates
print("Classifying HTT candidates...")
htt_rows = []
non_htt_rows = []
skipped = {"no_ks": 0, "no_qseqid": 0, "no_meta": 0, "no_mrca": 0, "no_threshold": 0}

for idx, r in te_ks.iterrows():
    if idx % 1000000 == 0 and idx > 0:
        print(f"  Processed {idx:,} pairs...")
    
    ks = r["ks"]
    if not np.isfinite(ks) or ks <= 0:
        skipped["no_ks"] += 1
        continue

    qseqid = str(r["qseqid"]) if pd.notna(r["qseqid"]) else ""
    sseqid = str(r["sseqid"]) if pd.notna(r["sseqid"]) else ""
    
    if not qseqid or qseqid == "nan" or not sseqid or sseqid == "nan":
        skipped["no_qseqid"] += 1
        continue

    key = qseqid + "||" + sseqid
    meta = pair_meta.get(key) or pair_meta.get(sseqid + "||" + qseqid)
    if not meta:
        skipped["no_meta"] += 1
        continue

    qcl = clade_map.get(meta["qgenome"], meta["qgenome"])
    scl = clade_map.get(meta["sgenome"], meta["sgenome"])

    mrca = find_mrca(qcl, scl)
    if not mrca:
        skipped["no_mrca"] += 1
        continue

    thr = busco_map.get(mrca)
    if thr is None or not np.isfinite(thr):
        skipped["no_threshold"] += 1
        continue

    is_htt = ks < thr
    
    rec = {
        "pair_id": r["pair_id"],
        "qseqid": qseqid,
        "sseqid": sseqid,
        "qgenome": meta["qgenome"],
        "sgenome": meta["sgenome"],
        "qclade": qcl,
        "sclade": scl,
        "mrca_node": mrca,
        "busco_ks_threshold": float(thr),
        "te_ks": float(ks),
        "ka": float(r["ka"]) if pd.notna(r["ka"]) and np.isfinite(r["ka"]) else None,
        "is_htt": is_htt
    }

    if is_htt:
        htt_rows.append(rec)
    else:
        non_htt_rows.append(rec)

print(f"")
print(f"HTT candidates: {len(htt_rows)}")
print(f"Non-HTT pairs: {len(non_htt_rows)}")
print(f"Skipped: {skipped}")

# Save results
if htt_rows:
    htt_df = pd.DataFrame(htt_rows)
    htt_df.to_csv(htt_dir / "htt_candidates.tsv", sep="\t", index=False)
    htt_df.groupby(["qclade", "sclade"]).size().reset_index(name="htt_count") \
        .to_csv(htt_dir / "htt_by_clade_pair.tsv", sep="\t", index=False)
    htt_df.groupby("mrca_node").size().reset_index(name="htt_count") \
        .to_csv(htt_dir / "htt_by_node.tsv", sep="\t", index=False)
    print(f"Wrote HTT candidates to {htt_dir / 'htt_candidates.tsv'}")

if non_htt_rows:
    non_df = pd.DataFrame(non_htt_rows)
    non_df.to_csv(htt_dir / "non_htt_pairs.tsv", sep="\t", index=False)
    print(f"Wrote non-HTT pairs to {htt_dir / 'non_htt_pairs.tsv'}")

# Summary
total = len(htt_rows) + len(non_htt_rows)
summary = [
    "HTT Detection Summary",
    "=" * 50,
    "",
    f"Total TE Ks pairs: {len(te_ks)}",
    f"Pairs with valid classification: {total}",
    f"HTT candidates: {len(htt_rows)}",
    f"Non-HTT pairs: {len(non_htt_rows)}",
]
if total > 0:
    summary.append(f"HTT percentage: {100 * len(htt_rows) / total:.4f}%")
if htt_rows:
    htt_ks = [r["te_ks"] for r in htt_rows]
    summary.append("")
    summary.append(f"HTT Ks range: {min(htt_ks):.4f} - {max(htt_ks):.4f}")
    summary.append(f"HTT Ks median: {np.median(htt_ks):.4f}")

summary.append("")
summary.append("Skipped breakdown:")
for k, v in skipped.items():
    summary.append(f"  {k}: {v}")

(summary_dir / "htt_summary.txt").write_text("\n".join(summary) + "\n")
print("Summary written.")
PYEOF

conda deactivate

echo ""
echo "STEP 7 complete at $(date)"

###############################################################################
### STEP 8: Generate plots
###############################################################################
echo ""
echo "### STEP 8: Generating plots ###"
echo ""

conda activate HTT_r

Rscript --vanilla - << 'REOF'
library(ggplot2)
library(dplyr)

res <- "/storage/zimmermf/HTT/results/06_results_htt_candidates"
ks_summary <- "/storage/zimmermf/HTT/results/03_results_ks_divergence/ks_calculations/node_ks_summary.tsv"

teks <- read.table(file.path(res, "ks_calculations", "te_ks_values.tsv"),
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("Loaded", nrow(teks), "TE Ks values\n")

httf <- file.path(res, "htt_candidates", "htt_candidates.tsv")
if (file.exists(httf) && file.size(httf) > 0) {
    htt <- read.table(httf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    teks$is_htt <- teks$pair_id %in% htt$pair_id
    cat("HTT candidates found:", sum(teks$is_htt), "\n")
} else {
    teks$is_htt <- FALSE
    cat("No HTT candidates file found\n")
}

teks <- teks %>% filter(is.finite(ks) & ks > 0)
cat("Valid Ks values:", nrow(teks), "\n")

busco <- read.table(ks_summary, header = TRUE, sep = "\t")
thr <- min(busco$ks_0.5_quantile, na.rm = TRUE)
cat("Min BUSCO threshold:", thr, "\n")

# Plot 1: Ks distribution
p <- ggplot(teks, aes(ks, fill = is_htt)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = thr, linetype = "dashed", color = "red") +
    scale_x_log10() +
    scale_fill_manual(
        values = c("FALSE" = "gray50", "TRUE" = "steelblue"),
        labels = c("Non-HTT", "HTT candidate")
    ) +
    labs(title = "TE Ks distribution", x = "Ks (log scale)", y = "Count", fill = "") +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave(file.path(res, "summary", "htt_ks_distribution.pdf"), p, width = 10, height = 6)
ggsave(file.path(res, "summary", "htt_ks_distribution.png"), p, width = 10, height = 6, dpi = 300)
cat("Saved Ks distribution plot\n")

# Plot 2: Heatmap if HTT candidates exist
cpf <- file.path(res, "htt_candidates", "htt_by_clade_pair.tsv")
if (file.exists(cpf) && file.size(cpf) > 0) {
    cp <- read.table(cpf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(cp) > 0) {
        p2 <- ggplot(cp, aes(x = qclade, y = sclade, fill = htt_count)) +
            geom_tile() +
            scale_fill_gradient(low = "white", high = "steelblue") +
            labs(title = "HTT candidates by clade pair", x = "Clade 1", y = "Clade 2", fill = "Count") +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                axis.text.y = element_text(size = 8)
            )
        ggsave(file.path(res, "summary", "htt_heatmap.pdf"), p2, width = 12, height = 10)
        ggsave(file.path(res, "summary", "htt_heatmap.png"), p2, width = 12, height = 10, dpi = 300)
        cat("Saved heatmap\n")
    }
}

# Plot 3: HTT by node
nf <- file.path(res, "htt_candidates", "htt_by_node.tsv")
if (file.exists(nf) && file.size(nf) > 0) {
    nd <- read.table(nf, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(nd) > 0) {
        p3 <- ggplot(nd, aes(x = reorder(mrca_node, htt_count), y = htt_count)) +
            geom_col(fill = "steelblue") +
            coord_flip() +
            labs(title = "HTT candidates by MRCA node", x = "MRCA Node", y = "HTT Count") +
            theme_minimal()
        ggsave(file.path(res, "summary", "htt_by_node.pdf"), p3, width = 10, height = 8)
        ggsave(file.path(res, "summary", "htt_by_node.png"), p3, width = 10, height = 8, dpi = 300)
        cat("Saved HTT by node plot\n")
    }
}

cat("Plots written\n")
REOF

conda deactivate

###############################################################################
### FINISH
###############################################################################
END=$(date +%s)
ELAPSED=$((END - START_TIME))
echo ""
echo "=========================================="
echo "Pipeline complete!"
echo "Runtime: $((ELAPSED / 3600))h $(((ELAPSED % 3600) / 60))m $((ELAPSED % 60))s"
echo ""
echo "Output files:"
echo "  ${RESULTS}/htt_candidates/htt_candidates.tsv"
echo "  ${RESULTS}/htt_candidates/non_htt_pairs.tsv"
echo "  ${RESULTS}/summary/htt_summary.txt"
echo "  ${RESULTS}/summary/htt_ks_distribution.pdf"
echo ""
echo "Backup location: ${TE_KS_BACKUP}"
echo "If anything goes wrong, restore with:"
echo "  cp ${TE_KS_BACKUP} ${TE_KS_OUT}"
echo "=========================================="