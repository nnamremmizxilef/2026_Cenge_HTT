#!/bin/bash
#SBATCH --job-name=htt_clustering
#SBATCH --qos normal
#SBATCH --account node
#SBATCH --partition node
#SBATCH --mail-user=felix.zimmermann@wsl.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --output /dev/null
#SBATCH --error /dev/null

set -euo pipefail

ROOT=/storage/zimmermf/HTT

TE_CLASS_RESULTS=${ROOT}/results/05_results_te_classification
HTT_RESULTS=${ROOT}/results/06_results_htt_candidates

RESULTS=${ROOT}/results/07_results_htt_clustering
LOGS=${ROOT}/logs/07_logs_htt_clustering

mkdir -p "${RESULTS}"/{within_clade/{fastas,db,blast},graphs,communities,merged_communities,summary,tmp}
mkdir -p "${LOGS}"

JOBTAG="${SLURM_JOB_ID}_htt_clustering"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "### HTT clustering started: $(date)"
echo "ROOT=${ROOT}"
echo "RESULTS=${RESULTS}"

# Inputs
HTT_FILE=${HTT_RESULTS}/htt_candidates/htt_candidates.tsv
CHAINED_HITS=${HTT_RESULTS}/chained_hits/chained_hits.tsv

TE_FASTA=${TE_CLASS_RESULTS}/all_tes_combined.fasta
TE_INFO=${TE_CLASS_RESULTS}/all_tes_info.tsv
CLADE_ASSIGN=${ROOT}/data/clade_assignment.tsv

# Paper-like BLAST params (within-clade)
BLAST_DBSIZE=10000000000
BLAST_EVALUE=10
BLAST_MAX_TARGETS=600000
MIN_LEN=300

# community merge threshold
MERGE_DENSITY=0.05

# Sanity checks
for f in "${HTT_FILE}" "${CHAINED_HITS}" "${TE_FASTA}" "${TE_INFO}" "${CLADE_ASSIGN}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo "HTT candidates: ${HTT_FILE}"
echo "Chained hits:   ${CHAINED_HITS}"
echo "TE fasta:       ${TE_FASTA}"

eval "$(conda shell.bash hook)"

###############################################################################
# STEP 1) Prepare per-clade TE FASTAs (only TEs that appear in HTT candidates)
###############################################################################
echo ""
echo "### STEP 1: Prepare per-clade TE FASTAs (HTT-only TEs) ###"
echo ""

conda activate HTT_python

python3 << 'PY'
import pandas as pd
from pathlib import Path
from collections import defaultdict

ROOT = Path("/storage/zimmermf/HTT")
RESULTS = ROOT / "results/07_results_htt_clustering"
HTT_FILE = ROOT / "results/06_results_htt_candidates/htt_candidates/htt_candidates.tsv"
TE_INFO  = ROOT / "results/05_results_te_classification/all_tes_info.tsv"
CLADE    = ROOT / "data/clade_assignment.tsv"
TE_FASTA = ROOT / "results/05_results_te_classification/all_tes_combined.fasta"

out_fasta_dir = RESULTS / "within_clade/fastas"
out_fasta_dir.mkdir(parents=True, exist_ok=True)

htt = pd.read_csv(HTT_FILE, sep="\t")
if htt.empty:
    raise SystemExit("No HTT candidates found; nothing to cluster.")

# TE -> genome
te_info = pd.read_csv(TE_INFO, sep="\t")
te_to_genome = dict(zip(te_info["te_id"], te_info["strain"]))

# genome(sample) -> clade
cl = pd.read_csv(CLADE, sep="\t")
cols = {c.lower(): c for c in cl.columns}
sample_col = cols.get("sample", cl.columns[0])
clade_col  = cols.get("clade",  cl.columns[1])
genome_to_clade = dict(zip(cl[sample_col], cl[clade_col]))

# gather all TE ids in HTT list
te_ids = set(htt["qseqid"]).union(set(htt["sseqid"]))

# clade -> set(te_ids)
clade_tes = defaultdict(set)
unknown = 0
for te in te_ids:
    g = te_to_genome.get(te)
    if g is None:
        unknown += 1
        continue
    c = genome_to_clade.get(g, g)
    clade_tes[c].add(te)

print(f"Unique HTT-involved TEs: {len(te_ids)}")
print(f"Clades with HTT-involved TEs: {len(clade_tes)}")
if unknown:
    print(f"WARNING: {unknown} TE ids not found in TE_INFO (will be ignored).")

# write a manifest
manifest = RESULTS / "summary/clade_te_manifest.tsv"
with open(manifest, "w") as out:
    out.write("clade\tn_tes\n")
    for c, s in sorted(clade_tes.items(), key=lambda x: (-len(x[1]), str(x[0]))):
        out.write(f"{c}\t{len(s)}\n")
print(f"Wrote manifest: {manifest}")

# build reverse map: te -> clade
te_to_clade = {}
for c, s in clade_tes.items():
    for te in s:
        te_to_clade[te] = c

need = set()
for s in clade_tes.values():
    need |= s

handles = {}
def get_handle(clade):
    if clade not in handles:
        fp = out_fasta_dir / f"{str(clade)}.fasta"
        handles[clade] = open(fp, "w")
    return handles[clade]

cur = None
buf = []
written = 0
with open(TE_FASTA) as f:
    for line in f:
        if line.startswith(">"):
            if cur is not None and cur in need:
                c = te_to_clade[cur]
                h = get_handle(c)
                h.write(f">{cur}\n{''.join(buf)}")
                written += 1
            cur = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if cur is not None and cur in need:
        c = te_to_clade[cur]
        h = get_handle(c)
        h.write(f">{cur}\n{''.join(buf)}")
        written += 1

for h in handles.values():
    h.close()

print(f"Wrote {written} sequences into per-clade FASTAs under {out_fasta_dir}")
PY

conda deactivate


###############################################################################
# STEP 2) Within-clade all-vs-all BLASTn (per clade)
###############################################################################
echo ""
echo "### STEP 2: Within-clade all-vs-all BLASTn ###"
echo ""

conda activate HTT_blast

FASTA_DIR="${RESULTS}/within_clade/fastas"
DB_DIR="${RESULTS}/within_clade/db"
BLAST_DIR="${RESULTS}/within_clade/blast"

mkdir -p "${DB_DIR}" "${BLAST_DIR}"

for fa in "${FASTA_DIR}"/*.fasta; do
  [ -e "$fa" ] || continue
  clade=$(basename "$fa" .fasta)

  nseq=$(grep -c "^>" "$fa" || true)
  if [ "${nseq}" -lt 2 ]; then
    echo "Skip clade ${clade}: only ${nseq} TE(s)"
    continue
  fi

  dbprefix="${DB_DIR}/${clade}"
  out="${BLAST_DIR}/${clade}_within.tsv"

  if [ -f "${out}" ]; then
    echo "Skip BLAST (exists): ${clade}"
    continue
  fi

  echo "BLAST within clade: ${clade} (nseq=${nseq})"

  makeblastdb -in "${fa}" -dbtype nucl -out "${dbprefix}" >/dev/null

  blastn \
    -query "${fa}" \
    -db "${dbprefix}" \
    -out "${out}" \
    -dbsize ${BLAST_DBSIZE} \
    -evalue ${BLAST_EVALUE} \
    -max_target_seqs ${BLAST_MAX_TARGETS} \
    -num_threads "${SLURM_CPUS_PER_TASK}" \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore evalue"

done

conda deactivate


###############################################################################
# STEP 2.5) Filter chained_hits to HTT-relevant pairs only (70GB -> few MB)
###############################################################################
echo ""
echo "### STEP 2.5: Filter chained_hits to HTT pairs only ###"
echo ""

CHAINED_FILTERED="${RESULTS}/tmp/chained_hits_filtered.tsv"

if [ -f "${CHAINED_FILTERED}" ]; then
  echo "Filtered chained_hits already exists, skipping..."
else
  # Extract unique TE IDs from HTT candidates (columns 2 and 3 = qseqid, sseqid)
  cut -f2,3 "${HTT_FILE}" | tail -n+2 | tr '\t' '\n' | sort -u > "${RESULTS}/tmp/htt_te_ids.txt"
  
  n_te_ids=$(wc -l < "${RESULTS}/tmp/htt_te_ids.txt")
  echo "Unique HTT TE IDs: ${n_te_ids}"
  
  # Filter chained_hits to only rows where both qseqid and sseqid are in HTT TEs
  echo "Filtering chained_hits (this may take a few minutes)..."
  awk -F'\t' 'NR==FNR {ids[$1]; next} FNR==1 || ($1 in ids && $2 in ids)' \
      "${RESULTS}/tmp/htt_te_ids.txt" \
      "${CHAINED_HITS}" > "${CHAINED_FILTERED}"
  
  echo "Filtered chained hits: $(wc -l < "${CHAINED_FILTERED}") rows"
fi


###############################################################################
# STEP 3) Cluster HTT candidates (igraph fastgreedy + merge communities)
###############################################################################
echo ""
echo "### STEP 3: Graph clustering per clade pair ###"
echo ""

conda activate HTT_python

python3 - "${MERGE_DENSITY}" "${MIN_LEN}" << 'PY'
import sys
import pandas as pd
from pathlib import Path
from collections import defaultdict
import math

MERGE_DENSITY = float(sys.argv[1])
MIN_LEN = int(sys.argv[2])

ROOT = Path("/storage/zimmermf/HTT")
RES  = ROOT / "results/07_results_htt_clustering"
HTT_FILE = ROOT / "results/06_results_htt_candidates/htt_candidates/htt_candidates.tsv"
CHAINED  = RES / "tmp/chained_hits_filtered.tsv"  # Use filtered version
WITHIN_DIR = RES / "within_clade/blast"

# --- igraph import ---
try:
    import igraph as ig
except Exception as e:
    raise SystemExit(f"ERROR: python-igraph not available. Install with: conda install -c conda-forge python-igraph\n{e}")

# Load HTTs
htt = pd.read_csv(HTT_FILE, sep="\t")
if htt.empty:
    raise SystemExit("No HTT candidates; nothing to cluster.")

required_cols = {"pair_id","qseqid","sseqid","qclade","sclade","mrca_node"}
missing = required_cols - set(htt.columns)
if missing:
    raise SystemExit(f"ERROR: HTT table missing columns: {missing}")

# Load between-TE identity from filtered chained hits
print(f"Loading filtered chained hits from: {CHAINED}")
ch = pd.read_csv(CHAINED, sep="\t")
print(f"Loaded {len(ch)} chained hit rows")
ch = ch[["qseqid","sseqid","pident"]].copy()

def norm_pair(a,b):
    return (a,b) if a <= b else (b,a)

ch["a"], ch["b"] = zip(*ch.apply(lambda r: norm_pair(r["qseqid"], r["sseqid"]), axis=1))
between_pid = dict(zip(zip(ch["a"], ch["b"]), ch["pident"]))

def get_between(row):
    a,b = norm_pair(row["qseqid"], row["sseqid"])
    return between_pid.get((a,b), math.nan)

htt["between_pident"] = htt.apply(get_between, axis=1)
htt = htt[pd.notna(htt["between_pident"])].copy()

print(f"HTT candidates with chained pident: {len(htt)}")

# Load within-clade pident
within = {}
for fp in WITHIN_DIR.glob("*_within.tsv"):
    clade = fp.name.replace("_within.tsv","")
    d = {}
    with open(fp) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            q,s,pid,L = parts[0], parts[1], parts[2], parts[3]
            if q == s:
                continue
            L = int(L)
            if L < MIN_LEN:
                continue
            a,b = norm_pair(q,s)
            pid = float(pid)
            if (a,b) not in d or pid > d[(a,b)]:
                d[(a,b)] = pid
    within[clade] = d

print(f"Loaded within-clade BLAST dicts: {len(within)} clades")

# Group HTTs by unordered clade-pair
def norm_clade_pair(a,b):
    return (a,b) if str(a) <= str(b) else (b,a)

htt["cladeA"], htt["cladeB"] = zip(*htt.apply(lambda r: norm_clade_pair(r["qclade"], r["sclade"]), axis=1))
groups = htt.groupby(["cladeA","cladeB"], sort=False)

out_comm_dir = RES / "communities"
out_merged_dir = RES / "merged_communities"
out_graph_dir = RES / "graphs"
out_sum_dir = RES / "summary"
for d in [out_comm_dir, out_merged_dir, out_graph_dir, out_sum_dir]:
    d.mkdir(exist_ok=True)

cluster_rows = []
summary_rows = []

# Global community counter for unique IDs across all clade pairs
global_community_id = 0

def within_pid(clade, te1, te2):
    d = within.get(clade)
    if not d:
        return 0.0
    a,b = norm_pair(te1, te2)
    return float(d.get((a,b), 0.0))

def build_edges(df, clA, clB):
    n = len(df)
    qmap = defaultdict(list)
    smap = defaultdict(list)
    for i,r in enumerate(df.itertuples(index=False)):
        qmap[r.qseqid].append(i)
        smap[r.sseqid].append(i)

    edges = set()

    # shared TE copy edges
    for lst in qmap.values():
        if len(lst) > 1:
            for i in range(len(lst)):
                for j in range(i+1, len(lst)):
                    edges.add((lst[i], lst[j]))
    for lst in smap.values():
        if len(lst) > 1:
            for i in range(len(lst)):
                for j in range(i+1, len(lst)):
                    edges.add((lst[i], lst[j]))

    # identity rule edges
    for i in range(n):
        qi, si, bi = df.iloc[i][["qseqid","sseqid","between_pident"]]
        for j in range(i+1, n):
            qj, sj, bj = df.iloc[j][["qseqid","sseqid","between_pident"]]
            wq = within_pid(clA, qi, qj)
            ws = within_pid(clB, si, sj)
            if max(float(bi), float(bj)) > max(wq, ws):
                edges.add((i,j))

    return list(edges)

def merge_communities(g, membership, density_thr=0.05):
    comm = defaultdict(list)
    for v, c in enumerate(membership):
        comm[c].append(v)

    edges = g.get_edgelist()
    inter = defaultdict(int)
    for a,b in edges:
        ca, cb = membership[a], membership[b]
        if ca == cb:
            continue
        if ca > cb:
            ca, cb = cb, ca
        inter[(ca,cb)] += 1

    parent = {c:c for c in comm.keys()}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a,b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for (ca,cb), ecount in inter.items():
        na, nb = len(comm[ca]), len(comm[cb])
        possible = na * nb
        dens = (ecount / possible) if possible else 0.0
        if dens > density_thr:
            union(ca,cb)

    remap = {}
    new_id = 0
    new_membership = []
    for v, c in enumerate(membership):
        rc = find(c)
        if rc not in remap:
            remap[rc] = new_id
            new_id += 1
        new_membership.append(remap[rc])

    return new_membership

for (clA, clB), df in groups:
    df = df.reset_index(drop=True)
    n_candidates = len(df)
    
    if n_candidates < 2:
        # singletons
        for r in df.itertuples(index=False):
            cluster_rows.append({
                "community_id": f"comm_{global_community_id:08d}",
                "pair_id": r.pair_id,
                "qseqid": r.qseqid,
                "sseqid": r.sseqid,
                "qclade": r.qclade,
                "sclade": r.sclade,
                "cladeA": clA,
                "cladeB": clB
            })
            global_community_id += 1
        summary_rows.append([clA, clB, n_candidates, 1, 1, 1, "singleton"])
        continue

    print(f"\nClustering clade pair: {clA} vs {clB} (n={n_candidates})")

    edges = build_edges(df, clA, clB)
    g = ig.Graph(n=n_candidates, edges=edges, directed=False)
    g.vs["pair_id"] = list(df["pair_id"])

    # save graph
    graph_path = out_graph_dir / f"{clA}_vs_{clB}.graphml"
    try:
        g.write_graphml(str(graph_path))
    except Exception:
        pass

    print(f"  Using fastgreedy algorithm")
    vc = g.community_fastgreedy()
    comm = vc.as_clustering()
    membership = list(comm.membership)
    algo_used = "fastgreedy"

    # write raw communities
    raw_path = out_comm_dir / f"{clA}_vs_{clB}_communities.tsv"
    pd.DataFrame({"pair_id": df["pair_id"], "community": membership}).to_csv(raw_path, sep="\t", index=False)

    # merge communities by density threshold
    merged = merge_communities(g, membership, density_thr=MERGE_DENSITY)

    merged_path = out_merged_dir / f"{clA}_vs_{clB}_clusters.tsv"
    pd.DataFrame({"pair_id": df["pair_id"], "cluster": merged}).to_csv(merged_path, sep="\t", index=False)

    # Map local cluster IDs to global community IDs
    local_to_global = {}
    for local_id in set(merged):
        local_to_global[local_id] = f"comm_{global_community_id:08d}"
        global_community_id += 1

    # Collect rows with full info for step 08
    for i, (pid, local_cid) in enumerate(zip(df["pair_id"], merged)):
        row = df.iloc[i]
        cluster_rows.append({
            "community_id": local_to_global[local_cid],
            "pair_id": pid,
            "qseqid": row["qseqid"],
            "sseqid": row["sseqid"],
            "qclade": row["qclade"],
            "sclade": row["sclade"],
            "cladeA": clA,
            "cladeB": clB
        })

    # summary
    sizes = pd.Series(merged).value_counts().sort_values(ascending=False)
    n_clusters = len(sizes)
    summary_rows.append([clA, clB, n_candidates, n_clusters, int(sizes.min()), int(sizes.max()), algo_used])

# Write combined output for step 08
clusters_df = pd.DataFrame(cluster_rows)
clusters_out = RES / "summary/htt_pairs_with_community.tsv"
clusters_df.to_csv(clusters_out, sep="\t", index=False)

# Also write simpler format
simple_out = RES / "summary/htt_clusters_all.tsv"
clusters_df[["cladeA","cladeB","pair_id","community_id"]].to_csv(simple_out, sep="\t", index=False)

summary_out = RES / "summary/htt_cluster_summary.tsv"
pd.DataFrame(summary_rows, columns=["cladeA","cladeB","n_candidates","n_clusters","min_cluster_size","max_cluster_size","algorithm"]).to_csv(summary_out, sep="\t", index=False)

print(f"\nWrote: {clusters_out}")
print(f"Wrote: {simple_out}")
print(f"Wrote: {summary_out}")
print(f"\nTotal communities: {global_community_id}")
PY

conda deactivate

echo ""
echo "### HTT clustering finished: $(date)"
echo "Results in: ${RESULTS}"
echo "Key outputs:"
echo "  - ${RESULTS}/summary/htt_pairs_with_community.tsv  <- for step 08"
echo "  - ${RESULTS}/summary/htt_clusters_all.tsv"
echo "  - ${RESULTS}/summary/htt_cluster_summary.tsv"