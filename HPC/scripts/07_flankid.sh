#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
GENOMES_DIR=${ROOT}/data/reference_genomes
HTT_RESULTS=${ROOT}/results/06_results_htt_candidates
RESULTS=${ROOT}/results/07_results_flankid
LOGS=${ROOT}/logs/07_logs_flankid

### input files from previous pipeline steps
IN_HTT=${HTT_RESULTS}/htt_candidates/htt_candidates.tsv
IN_NON=${HTT_RESULTS}/htt_candidates/non_htt_pairs.tsv

### flanking region parameters
# FLANK: size of upstream/downstream flanking regions to extract (bp)
FLANK=10000

### non-HTT matching mode
# per_family_equal: match equal number of non-HTT pairs per TE family
NON_MATCH_MODE="per_family_equal"

### TE alignment parameters (sensitive)
TE_NUCMER_MINMATCH=20
TE_NUCMER_MINCLUSTER=65

### flank alignment parameters (synteny-driven)
FLANK_NUCMER_MODE="--mum"
FLANK_NUCMER_MINMATCH=40
FLANK_NUCMER_MINCLUSTER=100
FLANK_MIN_ALN_LEN=200
TE_EXCLUDE_PAD=100



### create results/log folder(s)
mkdir -p "${RESULTS}"/{tmp,seqs,nucmer,summary}
mkdir -p "${LOGS}"



### setup logging
JOBTAG="${SLURM_JOB_ID:-local}_flankid"
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
echo "HTT candidates: ${IN_HTT}"
echo "Non-HTT pairs: ${IN_NON}"
echo "Genomes directory: ${GENOMES_DIR}"
echo "Results directory: ${RESULTS}"
echo ""
echo "Flanking parameters:"
echo "  Flank size: ${FLANK} bp"
echo "  TE exclude padding: ${TE_EXCLUDE_PAD} bp"
echo ""
echo "TE NUCmer parameters:"
echo "  Min match: ${TE_NUCMER_MINMATCH}"
echo "  Min cluster: ${TE_NUCMER_MINCLUSTER}"
echo ""
echo "Flank NUCmer parameters:"
echo "  Mode: ${FLANK_NUCMER_MODE}"
echo "  Min match: ${FLANK_NUCMER_MINMATCH}"
echo "  Min cluster: ${FLANK_NUCMER_MINCLUSTER}"
echo "  Min alignment length: ${FLANK_MIN_ALN_LEN}"



### sanity checks
for f in "${IN_HTT}" "${IN_NON}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        echo "Run 06_htt_candidates.sh first"
        exit 1
    fi
done

if [ ! -d "${GENOMES_DIR}" ]; then
    echo "ERROR: Genomes directory not found: ${GENOMES_DIR}"
    exit 1
fi

HTT_COUNT=$(tail -n +2 "${IN_HTT}" | wc -l)
NON_COUNT=$(tail -n +2 "${IN_NON}" | wc -l)
echo "Found ${HTT_COUNT} HTT candidates"
echo "Found ${NON_COUNT} non-HTT pairs"



### ===========================================================================
### STEP 1: Build matched pair list (HTT + matched non-HTT)
### ===========================================================================
echo ""
echo "### STEP 1: Building matched pair list ###"
echo ""

conda deactivate
conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"

PAIRLIST="${RESULTS}/summary/flankid_pairlist.tsv"

export IN_HTT IN_NON PAIRLIST

python3 << 'EOF'
import os
from pathlib import Path
import pandas as pd
import random

in_htt = os.environ.get("IN_HTT", "")
in_non = os.environ.get("IN_NON", "")
pairlist_out = os.environ.get("PAIRLIST", "")

htt = pd.read_csv(in_htt, sep="\t")
non = pd.read_csv(in_non, sep="\t")

need = {"qseqid", "sseqid", "qgenome", "sgenome", "is_htt"}
for name, df in [("HTT", htt), ("NON", non)]:
    miss = need - set(df.columns)
    if miss:
        raise SystemExit(f"ERROR: {name} missing columns: {sorted(miss)}")

def te_family(te_id: str) -> str:
    s = str(te_id)
    if "|" in s:
        s = s.split("|", 1)[1]
    if "::" in s:
        return s.split("::", 1)[0]
    return "UNKNOWN"

htt["family"] = htt["qseqid"].map(te_family)
non["family"] = non["qseqid"].map(te_family)

htt = htt.dropna(subset=["qseqid", "sseqid", "qgenome", "sgenome"]).copy()
non = non.dropna(subset=["qseqid", "sseqid", "qgenome", "sgenome"]).copy()

random.seed(1)
selected_non = []
target_counts = htt["family"].value_counts().to_dict()

for fam, k in target_counts.items():
    pool = non[non["family"] == fam]
    if len(pool) == 0:
        continue
    if len(pool) <= k:
        selected_non.append(pool)
    else:
        selected_non.append(pool.sample(n=k, replace=False, random_state=1))

non_sel = pd.concat(selected_non, ignore_index=True) if selected_non else non.iloc[0:0].copy()

htt["label"] = "HTT"
non_sel["label"] = "nonHTT"

pairs = pd.concat([htt, non_sel], ignore_index=True)
pairs = pairs[["label", "family", "qseqid", "sseqid", "qgenome", "sgenome"]].copy()
pairs.to_csv(pairlist_out, sep="\t", index=False)

print(f"HTT pairs: {len(htt)}")
print(f"Matched non-HTT pairs: {len(non_sel)}")
print(f"Total pairs: {len(pairs)}")
print(f"Wrote: {pairlist_out}")
EOF

echo ""
echo "STEP 1 complete."

conda deactivate



### ===========================================================================
### STEP 2: Extract TE+flanks sequences from genomes
### ===========================================================================
echo ""
echo "### STEP 2: Extracting TE+flanks sequences ###"
echo ""

conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"

export GENOMES_DIR RESULTS FLANK PAIRLIST

python3 << 'EOF'
import os
from pathlib import Path
import pandas as pd
import re
from collections import Counter
from functools import lru_cache

pairlist = Path(os.environ.get("PAIRLIST", ""))
genomes_dir = Path(os.environ.get("GENOMES_DIR", ""))
results_dir = Path(os.environ.get("RESULTS", ""))
out_seqs = results_dir / "seqs"
out_seqs.mkdir(parents=True, exist_ok=True)

FLANK = int(os.environ.get("FLANK", "10000"))

pat = re.compile(r"^(?P<genome>[^|]+)\|(?P<fam>[^:]+)::(?P<contig>[^:]+):(?P<start>\d+)-(?P<end>\d+)")

@lru_cache(maxsize=4)
def load_fasta_indexed(path: str):
    """Load FASTA with LRU cache to limit memory."""
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf).upper()
    return seqs

def get_genome_path(genome_name: str) -> Path:
    for ext in [".fasta", ".fa"]:
        p = genomes_dir / f"{genome_name}{ext}"
        if p.exists():
            return p
    for ext in [".fasta", ".fa"]:
        matches = [f for f in genomes_dir.glob(f"{genome_name}*{ext}") 
                   if not any(x in f.name for x in ['.bak', '.dict', '.fai', '.prep', '.ndb', '.nhr', '.nin'])]
        if matches:
            return matches[0]
    raise FileNotFoundError(f"Genome FASTA not found for {genome_name}")

def get_seq(genome_name: str, contig: str, start: int, end: int) -> str:
    fa_path = get_genome_path(genome_name)
    g = load_fasta_indexed(str(fa_path))
    
    if contig not in g:
        raise KeyError(f"Contig {contig} not found in {genome_name}")
    
    seq = g[contig]
    L = len(seq)
    s = max(1, start)
    e = min(L, end)
    return seq[s-1:e]

pairs = pd.read_csv(pairlist, sep="\t")
rows = []
fails = Counter()

for i, r in pairs.iterrows():
    qid = r["qseqid"]
    sid = r["sseqid"]
    mq = pat.match(str(qid))
    ms = pat.match(str(sid))
    
    if not mq or not ms:
        fails["bad_id_format"] += 1
        continue

    q = mq.groupdict()
    s = ms.groupdict()
    q0, q1 = int(q["start"]), int(q["end"])
    s0, s1 = int(s["start"]), int(s["end"])

    q_win0 = q0 - FLANK
    q_win1 = q1 + FLANK
    s_win0 = s0 - FLANK
    s_win1 = s1 + FLANK

    try:
        qseq = get_seq(q["genome"], q["contig"], q_win0, q_win1)
        sseq = get_seq(s["genome"], s["contig"], s_win0, s_win1)
    except FileNotFoundError:
        fails["missing_fasta"] += 1
        continue
    except KeyError:
        fails["missing_contig"] += 1
        continue
    except Exception:
        fails["other_extract_error"] += 1
        continue

    pair_id = f"pair_{i:09d}"
    qfa = out_seqs / f"{pair_id}_q.fa"
    sfa = out_seqs / f"{pair_id}_s.fa"

    # TE region within extracted window (1-based)
    q_te_start = max(1, FLANK + 1 + max(0, 1 - q_win0))
    q_te_end = q_te_start + (q1 - q0)
    s_te_start = max(1, FLANK + 1 + max(0, 1 - s_win0))
    s_te_end = s_te_start + (s1 - s0)

    with qfa.open("w") as out:
        out.write(f">{pair_id}|Q|{qid}|win:{q_win0}-{q_win1}|TE:{q_te_start}-{q_te_end}\n{qseq}\n")
    with sfa.open("w") as out:
        out.write(f">{pair_id}|S|{sid}|win:{s_win0}-{s_win1}|TE:{s_te_start}-{s_te_end}\n{sseq}\n")

    rows.append({
        "pair_id": pair_id,
        "label": r["label"],
        "family": r["family"],
        "qseqid": qid,
        "sseqid": sid,
        "qfa": str(qfa),
        "sfa": str(sfa),
        "q_te_start": q_te_start,
        "q_te_end": q_te_end,
        "s_te_start": s_te_start,
        "s_te_end": s_te_end,
        "q_len": len(qseq),
        "s_len": len(sseq),
    })

meta = pd.DataFrame(rows)
meta_out = results_dir / "summary" / "pair_meta.tsv"
meta.to_csv(meta_out, sep="\t", index=False)

print(f"Pairs total: {len(pairs)}")
print(f"Pairs with extracted windows: {len(meta)}")
if fails:
    print(f"Dropped pairs (reasons): {dict(fails)}")
print(f"Wrote: {meta_out}")
EOF

echo ""
echo "STEP 2 complete."

conda deactivate



### ===========================================================================
### STEP 3: Run NUCmer and compute coverage metrics
### ===========================================================================
echo ""
echo "### STEP 3: Running NUCmer and computing TE/flank metrics ###"
echo ""

conda activate HTT_mummer

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "nucmer version: $(nucmer --version 2>&1 | head -1)"

export RESULTS TE_NUCMER_MINMATCH TE_NUCMER_MINCLUSTER 
export FLANK_NUCMER_MODE FLANK_NUCMER_MINMATCH FLANK_NUCMER_MINCLUSTER 
export FLANK_MIN_ALN_LEN TE_EXCLUDE_PAD

python3 << 'EOF'
import os
from pathlib import Path
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

results_dir = Path(os.environ.get("RESULTS", ""))
meta_file = results_dir / "summary" / "pair_meta.tsv"
nucdir = results_dir / "nucmer"
nucdir.mkdir(exist_ok=True)

TE_minmatch = int(os.environ.get("TE_NUCMER_MINMATCH", "20"))
TE_mincluster = int(os.environ.get("TE_NUCMER_MINCLUSTER", "65"))
FL_mode = os.environ.get("FLANK_NUCMER_MODE", "--mum").strip()
FL_minmatch = int(os.environ.get("FLANK_NUCMER_MINMATCH", "40"))
FL_mincluster = int(os.environ.get("FLANK_NUCMER_MINCLUSTER", "100"))
FL_min_len = int(os.environ.get("FLANK_MIN_ALN_LEN", "200"))
TE_pad = int(os.environ.get("TE_EXCLUDE_PAD", "100"))
NCPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))

meta = pd.read_csv(meta_file, sep="\t")
if len(meta) == 0:
    raise SystemExit("No extracted pairs available; nothing to do.")

def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def union_covered(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    cov = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            cov += (cur_e - cur_s + 1)
            cur_s, cur_e = s, e
    cov += (cur_e - cur_s + 1)
    return cov

def intersect(intervals, a0, a1):
    out = []
    for s, e in intervals:
        s2 = max(s, a0)
        e2 = min(e, a1)
        if s2 <= e2:
            out.append((s2, e2))
    return out

def parse_coords(coords_path: Path):
    """Parse show-coords -rclT output."""
    alns = []
    with coords_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("[") or line.startswith("=") or "NUCMER" in line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                s1, e1 = int(parts[0]), int(parts[1])
                s2, e2 = int(parts[2]), int(parts[3])
                pid = float(parts[6])
            except ValueError:
                continue
            q0, q1 = min(s1, e1), max(s1, e1)
            t0, t1 = min(s2, e2), max(s2, e2)
            ori = 1 if (e1 - s1) * (e2 - s2) >= 0 else -1
            alns.append((q0, q1, t0, t1, pid, (q1 - q0 + 1), ori))
    return alns

def longest_collinear_chain_cov(alns, q_region, s_region, min_len=200):
    """Compute longest collinear chain coverage using DP."""
    q0r, q1r = q_region
    s0r, s1r = s_region

    keep = []
    for q0, q1, s0, s1, pid, lq, ori in alns:
        if lq < min_len:
            continue
        if q1 < q0r or q0 > q1r:
            continue
        if s1 < s0r or s0 > s1r:
            continue
        keep.append((q0, q1, s0, s1, ori))

    if not keep:
        return 0

    tot = {+1: 0, -1: 0}
    for q0, q1, s0, s1, ori in keep:
        tot[ori] += (q1 - q0 + 1)
    best_ori = +1 if tot[+1] >= tot[-1] else -1
    keep = [x for x in keep if x[4] == best_ori]
    if not keep:
        return 0

    keep.sort(key=lambda x: (x[0], x[1]))
    n = len(keep)
    
    dp = [0] * n
    prev = [-1] * n
    for i in range(n):
        q0, q1, s0, s1, ori = keep[i]
        dp[i] = (q1 - q0 + 1)
        for j in range(i):
            q0j, q1j, s0j, s1j, orj = keep[j]
            if s0 >= s0j and s1 >= s1j:
                cand = dp[j] + (q1 - q0 + 1)
                if cand > dp[i]:
                    dp[i] = cand
                    prev[i] = j

    end = max(range(n), key=lambda i: dp[i])
    chain = []
    while end != -1:
        chain.append(keep[end])
        end = prev[end]
    chain.reverse()

    q_intervals = [(max(q0, q0r), min(q1, q1r)) for q0, q1, _, _, _ in chain]
    q_intervals = [(a, b) for a, b in q_intervals if a <= b]
    return union_covered(q_intervals)

def process_pair(r):
    """Process a single pair: run nucmer and compute metrics."""
    pair_id = r["pair_id"]
    qfa = r["qfa"]
    sfa = r["sfa"]
    prefix_te = nucdir / f"{pair_id}.te"
    prefix_fl = nucdir / f"{pair_id}.fl"

    delta_te = str(prefix_te) + ".delta"
    coords_te = Path(str(prefix_te) + ".coords")
    delta_fl = str(prefix_fl) + ".delta"
    coords_fl = Path(str(prefix_fl) + ".coords")

    try:
        # TE: sensitive alignment
        run(["nucmer", "--maxmatch", "-l", str(TE_minmatch), "-c", str(TE_mincluster),
             "-p", str(prefix_te), qfa, sfa])
        with coords_te.open("w") as out:
            subprocess.run(["show-coords", "-rclT", delta_te], check=True, stdout=out, stderr=subprocess.DEVNULL)

        # FLANK: synteny-driven
        run(["nucmer", FL_mode, "-l", str(FL_minmatch), "-c", str(FL_mincluster),
             "-p", str(prefix_fl), qfa, sfa])
        with coords_fl.open("w") as out:
            subprocess.run(["show-coords", "-rclT", delta_fl], check=True, stdout=out, stderr=subprocess.DEVNULL)
    except Exception:
        return None

    q_len = int(r["q_len"])
    s_len = int(r["s_len"])
    q_te0, q_te1 = int(r["q_te_start"]), int(r["q_te_end"])
    s_te0, s_te1 = int(r["s_te_start"]), int(r["s_te_end"])

    te_alns = parse_coords(coords_te)
    if not te_alns:
        return None

    q_intervals = [(a, b) for a, b, _, _, _, _, _ in te_alns]
    s_intervals = [(c, d) for _, _, c, d, _, _, _ in te_alns]
    pids = [pid for _, _, _, _, pid, _, _ in te_alns]

    q_te_cov = union_covered(intersect(q_intervals, q_te0, q_te1)) / max(1, (q_te1 - q_te0 + 1))
    s_te_cov = union_covered(intersect(s_intervals, s_te0, s_te1)) / max(1, (s_te1 - s_te0 + 1))
    te_mean = (q_te_cov + s_te_cov) / 2.0

    fl_alns = parse_coords(coords_fl)

    def region_len(reg):
        a, b = reg
        return max(0, b - a + 1)

    q_up = (1, max(1, q_te0 - TE_pad - 1))
    q_dn = (min(q_len, q_te1 + TE_pad + 1), q_len)
    s_up = (1, max(1, s_te0 - TE_pad - 1))
    s_dn = (min(s_len, s_te1 + TE_pad + 1), s_len)

    q_up_cov_bp = longest_collinear_chain_cov(fl_alns, q_up, s_up, min_len=FL_min_len) if region_len(q_up) else 0
    q_dn_cov_bp = longest_collinear_chain_cov(fl_alns, q_dn, s_dn, min_len=FL_min_len) if region_len(q_dn) else 0
    
    fl_alns_swap = [(c, d, a, b, pid, lq, ori) for a, b, c, d, pid, lq, ori in fl_alns]
    s_up_cov_bp = longest_collinear_chain_cov(fl_alns_swap, s_up, q_up, min_len=FL_min_len) if region_len(s_up) else 0
    s_dn_cov_bp = longest_collinear_chain_cov(fl_alns_swap, s_dn, q_dn, min_len=FL_min_len) if region_len(s_dn) else 0

    q_up_cov = q_up_cov_bp / max(1, region_len(q_up)) if region_len(q_up) else 0.0
    q_dn_cov = q_dn_cov_bp / max(1, region_len(q_dn)) if region_len(q_dn) else 0.0
    s_up_cov = s_up_cov_bp / max(1, region_len(s_up)) if region_len(s_up) else 0.0
    s_dn_cov = s_dn_cov_bp / max(1, region_len(s_dn)) if region_len(s_dn) else 0.0

    q_flank_synt = (q_up_cov + q_dn_cov) / 2.0
    s_flank_synt = (s_up_cov + s_dn_cov) / 2.0
    flank_synt_mean = (q_flank_synt + s_flank_synt) / 2.0

    htt_score = te_mean - flank_synt_mean

    return {
        "pair_id": pair_id,
        "label": r["label"],
        "family": r["family"],
        "qseqid": r["qseqid"],
        "sseqid": r["sseqid"],
        "q_te_cov": q_te_cov,
        "s_te_cov": s_te_cov,
        "te_cov_mean": te_mean,
        "q_flank_synt_cov": q_flank_synt,
        "s_flank_synt_cov": s_flank_synt,
        "flank_synt_cov_mean": flank_synt_mean,
        "q_up_chain_bp": q_up_cov_bp,
        "q_dn_chain_bp": q_dn_cov_bp,
        "s_up_chain_bp": s_up_cov_bp,
        "s_dn_chain_bp": s_dn_cov_bp,
        "mean_pid": sum(pids) / len(pids) if pids else None,
        "htt_score_te_minus_flank": htt_score,
    }

# parallel processing
rows = []
print(f"Processing {len(meta)} pairs with {NCPUS} threads...")
with ThreadPoolExecutor(max_workers=NCPUS) as executor:
    futures = {executor.submit(process_pair, r): r["pair_id"] for _, r in meta.iterrows()}
    done = 0
    for future in as_completed(futures):
        done += 1
        if done % 500 == 0:
            print(f"  Processed {done}/{len(meta)} pairs...")
        result = future.result()
        if result is not None:
            rows.append(result)

out = pd.DataFrame(rows)
out_tsv = results_dir / "summary" / "flankid_results.tsv"
out.to_csv(out_tsv, sep="\t", index=False)
print(f"Computed results for pairs: {len(out)}")
print(f"Wrote: {out_tsv}")
EOF

echo ""
echo "STEP 3 complete."

conda deactivate



### ===========================================================================
### STEP 4: Generate summary stats
### ===========================================================================
echo ""
echo "### STEP 4: Generating summary stats ###"
echo ""

conda activate HTT_r

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "R version: $(R --version | head -1)"

export RESULTS

Rscript --vanilla << 'REOF'
library(dplyr)

results_dir <- Sys.getenv("RESULTS")
infile <- file.path(results_dir, "summary/flankid_results.tsv")
outdir <- file.path(results_dir, "summary")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df %>% filter(is.finite(te_cov_mean), is.finite(flank_synt_cov_mean))

# Summary stats
summary_stats <- df %>%
  group_by(label) %>%
  summarise(
    n = n(),
    te_cov_median = median(te_cov_mean),
    te_cov_mean = mean(te_cov_mean),
    flank_synt_median = median(flank_synt_cov_mean),
    flank_synt_mean = mean(flank_synt_cov_mean),
    htt_score_median = median(htt_score_te_minus_flank),
    htt_score_mean = mean(htt_score_te_minus_flank),
    .groups = "drop"
  )
write.table(summary_stats, file.path(outdir, "summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Summary written to:", outdir, "\n")
REOF

echo ""
echo "STEP 4 complete."

conda deactivate



### ===========================================================================
### STEP 5: Filter strict HTT candidates (TE identity >75%, flank synteny thresholds)
### ===========================================================================
echo ""
echo "### STEP 5: Filtering strict HTT candidates ###"
echo ""

conda activate HTT_python
echo "Running in environment: ${CONDA_DEFAULT_ENV}"

export RESULTS IN_HTT

python3 << 'EOF'
import os
from pathlib import Path
import pandas as pd

results_dir = Path(os.environ.get("RESULTS", ""))
in_htt = os.environ.get("IN_HTT", "")

infile = results_dir / "summary" / "flankid_results.tsv"
outdir = results_dir / "summary"
outdir.mkdir(parents=True, exist_ok=True)

if not infile.exists():
    raise SystemExit(f"ERROR: missing {infile}")

# Load flankid results and original HTT candidates
df_flank = pd.read_csv(infile, sep="\t")
df_orig = pd.read_csv(in_htt, sep="\t")

# Merge on qseqid+sseqid to get original pair_id and all columns
# Keep original columns (no suffix), mark flankid duplicates for removal
df = df_flank.merge(
    df_orig,
    on=["qseqid", "sseqid"],
    how="left",
    suffixes=("_flankid", "")
)

# Drop flankid-generated columns that duplicate originals
cols_to_drop = [c for c in df.columns if c.endswith("_flankid")]
df = df.drop(columns=cols_to_drop)

print(f"Merged flankid results with original HTT candidates")
print(f"  Flankid rows: {len(df_flank)}")
print(f"  Merged rows:  {len(df)}")

# Sequential filtering
n_all = len(df_flank)
df_htt = df[df["label"] == "HTT"].copy()
n_htt = len(df_htt)

df_htt["mean_pid"] = pd.to_numeric(df_htt["mean_pid"], errors="coerce")
df_te = df_htt[df_htt["mean_pid"] > 75.0].copy()
n_te = len(df_te)

df_te["flank_synt_cov_mean"] = pd.to_numeric(df_te["flank_synt_cov_mean"], errors="coerce")
eps = 1e-12

df_flank_0 = df_te[df_te["flank_synt_cov_mean"].fillna(1.0) <= eps].copy()
df_flank_1 = df_te[df_te["flank_synt_cov_mean"].fillna(1.0) <= 0.01].copy()
df_flank_5 = df_te[df_te["flank_synt_cov_mean"].fillna(1.0) <= 0.05].copy()

# Write filtered candidates (all columns preserved, original pair_id)
df_flank_0.to_csv(outdir / "htt_candidates_flank_0pct_teid_gt75.tsv", sep="\t", index=False)
df_flank_1.to_csv(outdir / "htt_candidates_flank_le1pct_teid_gt75.tsv", sep="\t", index=False)
df_flank_5.to_csv(outdir / "htt_candidates_flank_le5pct_teid_gt75.tsv", sep="\t", index=False)

# Build summary
summary_rows = [
    {"step": "all_pairs", "n_pairs": n_all, "filtered_out": 0},
    {"step": "HTT_only", "n_pairs": n_htt, "filtered_out": n_all - n_htt},
    {"step": "HTT_TE_identity_gt75", "n_pairs": n_te, "filtered_out": n_htt - n_te},
    {"step": "HTT_TEid_gt75_flank_0pct", "n_pairs": len(df_flank_0), "filtered_out": n_te - len(df_flank_0)},
    {"step": "HTT_TEid_gt75_flank_le1pct", "n_pairs": len(df_flank_1), "filtered_out": n_te - len(df_flank_1)},
    {"step": "HTT_TEid_gt75_flank_le5pct", "n_pairs": len(df_flank_5), "filtered_out": n_te - len(df_flank_5)},
]

summary = pd.DataFrame(summary_rows)
out_summary = outdir / "htt_candidate_filtering_summary.tsv"
summary.to_csv(out_summary, sep="\t", index=False)

print("\nFiltering summary:")
print(summary.to_string(index=False))
print(f"\nNote: flank thresholds are parallel filters from HTT_TE_identity_gt75")
print(f"\nWrote: {out_summary}")
EOF

echo ""
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
echo "  Pair list:                  ${RESULTS}/summary/flankid_pairlist.tsv"
echo "  Pair metadata:              ${RESULTS}/summary/pair_meta.tsv"
echo "  FlankID results:            ${RESULTS}/summary/flankid_results.tsv"
echo "  Summary statistics:         ${RESULTS}/summary/summary_stats.tsv"
echo "  Filtering summary:          ${RESULTS}/summary/htt_candidate_filtering_summary.tsv"
echo "  Strict HTT (0% flank):      ${RESULTS}/summary/htt_candidates_flank_0pct_teid_gt75.tsv"
echo "  Strict HTT (≤1% flank):     ${RESULTS}/summary/htt_candidates_flank_le1pct_teid_gt75.tsv"
echo "  Strict HTT (≤5% flank):     ${RESULTS}/summary/htt_candidates_flank_le5pct_teid_gt75.tsv"
echo "  Extracted sequences:        ${RESULTS}/seqs/"
echo ""