#!/bin/bash

set -euo pipefail

### paths
ROOT=/path/to/HTT
export ROOT

GENOMES_DIR=${ROOT}/data/reference_genomes

HTT_RESULTS=${ROOT}/results/06_results_htt_candidates
IN_HTT=${HTT_RESULTS}/htt_candidates/htt_candidates.tsv
IN_NON=${HTT_RESULTS}/htt_candidates/non_htt_pairs.tsv

RESULTS=${ROOT}/results/09_results_flankid
LOGS=${ROOT}/logs/09_logs_flankid

### params
FLANK=10000                # 10kb up/down
NON_MATCH_MODE="per_family_equal"

# TE alignment (sensitive)
TE_NUCMER_MINMATCH=20
TE_NUCMER_MINCLUSTER=65

# FLANK alignment (synteny-driven)
FLANK_NUCMER_MODE="--mum"
FLANK_NUCMER_MINMATCH=40
FLANK_NUCMER_MINCLUSTER=100
FLANK_MIN_ALN_LEN=200
TE_EXCLUDE_PAD=100

### dirs
mkdir -p "${RESULTS}"/{tmp,seqs,nucmer,summary}
mkdir -p "${LOGS}"

### logging
JOBTAG="${SLURM_JOB_ID:-local}_flankid"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "### 09_flankid start: $(date)"
echo "IN_HTT: ${IN_HTT}"
echo "IN_NON: ${IN_NON}"
echo "GENOMES_DIR: ${GENOMES_DIR}"
echo "FLANK: ${FLANK}"
echo "RESULTS: ${RESULTS}"

### sanity
for f in "${IN_HTT}" "${IN_NON}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing input: $f"
    exit 1
  fi
done
if [ ! -d "${GENOMES_DIR}" ]; then
  echo "ERROR: missing genomes dir: ${GENOMES_DIR}"
  exit 1
fi

eval "$(conda shell.bash hook)"
conda activate HTT_python

###############################################################################
# STEP 1: build pair list = all HTT + matched non-HTT (by TE family)
###############################################################################
echo ""
echo "### STEP 1: Build matched pair list (HTT + matched non-HTT) ###"
echo ""

PAIRLIST="${RESULTS}/summary/flankid_pairlist.tsv"

export IN_HTT IN_NON PAIRLIST

python3 << 'PY'
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
PY



###############################################################################
# STEP 2: Extract TE+flanks sequences from the two genomes
###############################################################################
echo ""
echo "### STEP 2: Extract TE+flanks sequences (per pair) ###"
echo ""

export GENOMES_DIR RESULTS FLANK PAIRLIST

python3 << 'PY'
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

# Fixed regex: removed $ anchor to handle trailing ::contig:start-end(strand) suffix
pat = re.compile(r"^(?P<genome>[^|]+)\|(?P<fam>[^:]+)::(?P<contig>[^:]+):(?P<start>\d+)-(?P<end>\d+)")

@lru_cache(maxsize=4)
def load_fasta_indexed(path: str):
    """Load FASTA with LRU cache to limit memory (keeps max 4 genomes)."""
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
    # Try exact match first
    for ext in [".fasta", ".fa"]:
        p = genomes_dir / f"{genome_name}{ext}"
        if p.exists():
            return p
    # Try glob pattern for longer names (e.g., Cenge1005 -> Cenge1005_1_AssemblyScaffolds_*.fasta)
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
    except FileNotFoundError as e:
        fails["missing_fasta"] += 1
        continue
    except KeyError as e:
        fails["missing_contig"] += 1
        continue
    except Exception as e:
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
PY



###############################################################################
# STEP 3: Run nucmer for each pair + compute coverage
###############################################################################
echo ""
echo "### STEP 3: NUCmer + compute TE / synteny-flank mapping rates ###"
echo ""

export RESULTS TE_NUCMER_MINMATCH TE_NUCMER_MINCLUSTER FLANK_NUCMER_MODE FLANK_NUCMER_MINMATCH FLANK_NUCMER_MINCLUSTER FLANK_MIN_ALN_LEN TE_EXCLUDE_PAD

python3 << 'PY'
import os
from pathlib import Path
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys

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
    """Parse show-coords -rclT output, skipping header lines robustly."""
    alns = []
    with coords_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Skip header lines
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

    # Choose best orientation
    tot = {+1: 0, -1: 0}
    for q0, q1, s0, s1, ori in keep:
        tot[ori] += (q1 - q0 + 1)
    best_ori = +1 if tot[+1] >= tot[-1] else -1
    keep = [x for x in keep if x[4] == best_ori]
    if not keep:
        return 0

    keep.sort(key=lambda x: (x[0], x[1]))
    n = len(keep)
    
    # DP for longest monotone chain
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

    # Traceback
    end = max(range(n), key=lambda i: dp[i])
    chain = []
    while end != -1:
        chain.append(keep[end])
        end = prev[end]
    chain.reverse()

    # Union coverage restricted to q_region
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
    except Exception as e:
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

    # Flank regions excluding TE ± pad
    q_up = (1, max(1, q_te0 - TE_pad - 1))
    q_dn = (min(q_len, q_te1 + TE_pad + 1), q_len)
    s_up = (1, max(1, s_te0 - TE_pad - 1))
    s_dn = (min(s_len, s_te1 + TE_pad + 1), s_len)

    q_up_cov_bp = longest_collinear_chain_cov(fl_alns, q_up, s_up, min_len=FL_min_len) if region_len(q_up) else 0
    q_dn_cov_bp = longest_collinear_chain_cov(fl_alns, q_dn, s_dn, min_len=FL_min_len) if region_len(q_dn) else 0
    
    # Swap q/s for subject-based calculation
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

# Parallel processing
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
PY

conda deactivate


###############################################################################
# STEP 4: Plots
###############################################################################
echo ""
echo "### STEP 4: Plot distributions (HTT vs non-HTT) ###"
echo ""

conda activate HTT_r

Rscript --vanilla - << 'RS'
library(ggplot2)
library(dplyr)

root <- Sys.getenv("ROOT", "/path/to/HTT")
infile <- file.path(root, "results/09_results_flankid/summary/flankid_results.tsv")
outdir <- file.path(root, "results/09_results_flankid/summary")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df <- df %>% filter(is.finite(te_cov_mean), is.finite(flank_synt_cov_mean))

# Color palette
cols <- c("HTT" = "#E64B35", "nonHTT" = "#4DBBD5")

# 1. Flank synteny histogram + density
p1 <- ggplot(df, aes(x = flank_synt_cov_mean, fill = label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.5, position = "identity") +
  geom_density(aes(color = label), linewidth = 0.8, fill = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = "Flank synteny coverage (collinear chain)", y = "Density",
       title = "Flanking region synteny: HTT vs non-HTT") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "flank_synt_cov_hist.png"), p1, width = 10, height = 6, dpi = 300)

# 2. TE coverage histogram + density
p2 <- ggplot(df, aes(x = te_cov_mean, fill = label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.5, position = "identity") +
  geom_density(aes(color = label), linewidth = 0.8, fill = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = "TE coverage mean", y = "Density", title = "TE mapping rate: HTT vs non-HTT") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "te_cov_hist.png"), p2, width = 10, height = 6, dpi = 300)

# 3. HTT score histogram + density
p3 <- ggplot(df, aes(x = htt_score_te_minus_flank, fill = label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.5, position = "identity") +
  geom_density(aes(color = label), linewidth = 0.8, fill = NA) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = "HTT score (TE_cov - flank_synt_cov)", y = "Density",
       title = "HTT score distribution") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "htt_score_hist.png"), p3, width = 10, height = 6, dpi = 300)

# 4. Scatter plot with marginal densities
p4 <- ggplot(df, aes(x = te_cov_mean, y = flank_synt_cov_mean, color = label)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = cols) +
  labs(x = "TE coverage mean", y = "Flank synteny coverage mean",
       title = "TE vs synteny-flank mapping") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "te_vs_flank_synt_scatter.png"), p4, width = 8, height = 7, dpi = 300)

# 5. Violin plot comparison
df_long <- df %>%
  tidyr::pivot_longer(cols = c(te_cov_mean, flank_synt_cov_mean),
                      names_to = "metric", values_to = "coverage") %>%
  mutate(metric = recode(metric,
                         "te_cov_mean" = "TE",
                         "flank_synt_cov_mean" = "Flank synteny"))

p5 <- ggplot(df_long, aes(x = label, y = coverage, fill = label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~metric) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Coverage", title = "Coverage distributions by category") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(outdir, "coverage_violin.png"), p5, width = 10, height = 6, dpi = 300)

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

cat("Plots and summary written to:", outdir, "\n")
RS

conda deactivate

echo ""
echo "### 09_flankid done: $(date)"
echo "Key outputs:"
echo "  ${RESULTS}/summary/flankid_pairlist.tsv"
echo "  ${RESULTS}/summary/pair_meta.tsv"
echo "  ${RESULTS}/summary/flankid_results.tsv"
echo "  ${RESULTS}/summary/summary_stats.tsv"
echo "  ${RESULTS}/summary/*.png"
