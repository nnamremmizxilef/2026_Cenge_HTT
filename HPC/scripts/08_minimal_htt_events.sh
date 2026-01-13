#!/bin/bash
#
#SBATCH --job-name=min_htt_events
#SBATCH --qos hprio
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

### paths
ROOT=/storage/zimmermf/HTT

TREE_RESULTS=${ROOT}/results/01_results_phylo_tree
KS_RESULTS=${ROOT}/results/03_results_ks_divergence
CLUST_RESULTS=${ROOT}/results/07_results_htt_clustering

RESULTS=${ROOT}/results/08_results_minimal_htt_events
LOGS=${ROOT}/logs/08_logs_minimal_htt_events

# inputs
TREE_FILE=${TREE_RESULTS}/phylogeny_rooted.treefile
NODES_INFO=${KS_RESULTS}/node_analyses/nodes_info.json
CLADE_ASSIGN=${ROOT}/data/clade_assignment.tsv

# Output from step 07 with columns: community_id, qseqid, sseqid, qclade, sclade, ...
CLUSTERED_HTT_TSV=${CLUST_RESULTS}/summary/htt_pairs_with_community.tsv


### create dirs
mkdir -p "${RESULTS}"/{network,clusters,events,trimmed_trees,summary,tmp}
mkdir -p "${LOGS}"

### logging
JOBTAG="${SLURM_JOB_ID}_minimal_htt_events"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "### 08_minimal_htt_events starting: $(date)"
echo "TREE_FILE: ${TREE_FILE}"
echo "NODES_INFO: ${NODES_INFO}"
echo "CLADE_ASSIGN: ${CLADE_ASSIGN}"
echo "CLUSTERED_HTT_TSV: ${CLUSTERED_HTT_TSV}"
echo "RESULTS: ${RESULTS}"

### sanity
for f in "${TREE_FILE}" "${NODES_INFO}" "${CLADE_ASSIGN}" "${CLUSTERED_HTT_TSV}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing input: $f"
    exit 1
  fi
done

eval "$(conda shell.bash hook)"
conda activate HTT_python_ete3


###############################################################################
# STEP 1: Build "community network" (nodes = communities, edges if share TE)
#         Single-linkage clustering = connected components
###############################################################################
echo ""
echo "### STEP 1: Build network of communities and connected components (single linkage) ###"
echo ""

python3 - "${CLUSTERED_HTT_TSV}" "${RESULTS}" << 'PY'
import sys
from pathlib import Path
import pandas as pd
import networkx as nx

clustered_tsv = Path(sys.argv[1])
results_dir = Path(sys.argv[2])
outdir = results_dir / "network"
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(clustered_tsv, sep="\t")
required = {"community_id","qseqid","sseqid","qclade","sclade"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"ERROR: {clustered_tsv} missing columns: {sorted(missing)}")

df["qseqid"] = df["qseqid"].astype(str)
df["sseqid"] = df["sseqid"].astype(str)
df["community_id"] = df["community_id"].astype(str)

# Build TE -> communities index
te_to_comms = {}
for tecol in ["qseqid","sseqid"]:
    for te, sub in df.groupby(tecol):
        te_to_comms.setdefault(te, set()).update(sub["community_id"].unique())

G = nx.Graph()
G.add_nodes_from(df["community_id"].unique())

for te, comms in te_to_comms.items():
    comms = list(comms)
    if len(comms) < 2:
        continue
    for i in range(len(comms)):
        for j in range(i+1, len(comms)):
            a, b = comms[i], comms[j]
            if a == b:
                continue
            if not G.has_edge(a, b):
                G.add_edge(a, b, shared_te=set())
            G.edges[a, b]["shared_te"].add(te)

print(f"Communities: {G.number_of_nodes():,}")
print(f"Edges (shared TE): {G.number_of_edges():,}")

# Single linkage clustering = connected components
components = list(nx.connected_components(G))
print(f"Connected components (single-linkage clusters): {len(components):,}")

# Write mapping: community -> component
rows = []
for cid, comp in enumerate(components, start=1):
    for comm in comp:
        rows.append((comm, f"component_{cid:06d}", len(comp)))

map_df = pd.DataFrame(rows, columns=["community_id","component_id","component_n_communities"])
map_df.to_csv(outdir / "community_to_component.tsv", sep="\t", index=False)

# Component summary
comp_sizes = map_df.groupby("component_id")["community_id"].count().reset_index(name="n_communities")
comp_sizes.to_csv(outdir / "component_sizes.tsv", sep="\t", index=False)

print("Wrote:")
print(" -", outdir / "community_to_component.tsv")
print(" -", outdir / "component_sizes.tsv")
PY


###############################################################################
# STEP 2: For each component, trim tree with ete3 and count minimal HTT events
###############################################################################
echo ""
echo "### STEP 2: Minimal HTT events per component (with ete3 tree trimming) ###"
echo ""

python3 - "${TREE_FILE}" "${NODES_INFO}" "${CLADE_ASSIGN}" "${CLUSTERED_HTT_TSV}" "${RESULTS}" << 'PY'
import sys
from pathlib import Path
import json
import pandas as pd
from collections import defaultdict

tree_file = Path(sys.argv[1])
nodes_info_file = Path(sys.argv[2])
clade_file = Path(sys.argv[3])
clustered_tsv = Path(sys.argv[4])
results_dir = Path(sys.argv[5])

map_file = results_dir / "network" / "community_to_component.tsv"
out_events = results_dir / "events"
out_trees = results_dir / "trimmed_trees"
out_summary = results_dir / "summary"
out_events.mkdir(parents=True, exist_ok=True)
out_trees.mkdir(parents=True, exist_ok=True)
out_summary.mkdir(parents=True, exist_ok=True)

# Try to import ete3
try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False
    print("WARNING: ete3 not available. Tree trimming will be skipped.")
    print("Install with: conda install -c etetoolkit ete3")

# --- Load clade map: genome -> clade
clade_map = {}
with open(clade_file) as f:
    header = f.readline().strip().split("\t")
    col = {h.lower():i for i,h in enumerate(header)}
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) < 2:
            continue
        sample = fields[col.get("sample", 0)]
        clade = fields[col.get("clade", 1)]
        clade_map[sample] = clade

# --- Load species tree
if HAS_ETE3:
    # Read tree and simplify tip labels
    tree_str = tree_file.read_text().strip()
    full_tree = Tree(tree_str, format=1)
    
    # Simplify tip names to match clade_map keys
    for leaf in full_tree.iter_leaves():
        # Remove _1_AssemblyScaffolds... suffix
        simple_name = leaf.name.split("_1_AssemblyScaffolds")[0]
        if simple_name.endswith("_1"):
            simple_name = simple_name[:-2]
        leaf.name = simple_name
    
    print(f"Loaded tree with {len(full_tree)} leaves")

# --- Load nodes_info
nodes_info = json.loads(nodes_info_file.read_text())

# Build node_to_clades sets
node_to_clades = {}
node_size = {}
for node in nodes_info:
    clades = set()
    for leaf in node["all_leaves"]:
        # Simplify leaf name
        simple_leaf = leaf.split("_1_AssemblyScaffolds")[0]
        if simple_leaf.endswith("_1"):
            simple_leaf = simple_leaf[:-2]
        clades.add(clade_map.get(simple_leaf, simple_leaf))
    node_to_clades[node["node_name"]] = clades
    node_size[node["node_name"]] = len(clades)

def find_mrca_node(c1, c2):
    """Find MRCA node containing both clades."""
    candidates = [n for n, cs in node_to_clades.items() if c1 in cs and c2 in cs]
    if not candidates:
        return None
    return min(candidates, key=lambda n: node_size[n])

# --- Load HTT candidates with community assignments
df = pd.read_csv(clustered_tsv, sep="\t")
required = {"community_id","qseqid","sseqid","qclade","sclade"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"ERROR: {clustered_tsv} missing columns: {sorted(missing)}")

df["qclade"] = df["qclade"].astype(str)
df["sclade"] = df["sclade"].astype(str)
df["community_id"] = df["community_id"].astype(str)

# Add component assignment
mp = pd.read_csv(map_file, sep="\t")
mp["community_id"] = mp["community_id"].astype(str)
df = df.merge(mp[["community_id","component_id"]], on="community_id", how="left")

if df["component_id"].isna().any():
    n_missing = df["component_id"].isna().sum()
    print(f"WARNING: {n_missing} rows without component assignment")
    df = df[df["component_id"].notna()].copy()

# Deduplicate to unique community entries per clade pair
df_comm = df.drop_duplicates(subset=["community_id","qclade","sclade"]).copy()

def norm_pair(a,b):
    return tuple(sorted((a,b)))

df_comm["clade_pair"] = [norm_pair(a,b) for a,b in zip(df_comm["qclade"], df_comm["sclade"])]

# --- Per component: trim tree and count minimal events
component_rows = []
node_event_counts = defaultdict(int)
total_components = df_comm["component_id"].nunique()

print(f"Processing {total_components} components...")

for comp_idx, (comp_id, sub) in enumerate(df_comm.groupby("component_id")):
    if (comp_idx + 1) % 100 == 0:
        print(f"  Processed {comp_idx + 1}/{total_components} components")
    
    # Get all clades in this component
    clades_in_component = set(sub["qclade"]).union(set(sub["sclade"]))
    
    # Get all genomes (samples) belonging to these clades
    samples_in_component = set()
    for sample, clade in clade_map.items():
        if clade in clades_in_component:
            samples_in_component.add(sample)
    
    # Also add clades that are their own samples (uncollapsed)
    for clade in clades_in_component:
        if clade not in [c for c in clade_map.values()]:
            # This clade might be a sample name itself
            samples_in_component.add(clade)
    
    # Trim tree to only these samples (if ete3 available)
    if HAS_ETE3 and len(samples_in_component) >= 2:
        try:
            # Find leaves that match our samples
            leaves_to_keep = []
            for leaf in full_tree.iter_leaves():
                if leaf.name in samples_in_component:
                    leaves_to_keep.append(leaf.name)
            
            if len(leaves_to_keep) >= 2:
                trimmed = full_tree.copy()
                trimmed.prune(leaves_to_keep, preserve_branch_length=True)
                
                # Save trimmed tree (optional, for large datasets might skip)
                if total_components <= 1000:  # Only save for manageable number
                    tree_path = out_trees / f"{comp_id}.nwk"
                    trimmed.write(outfile=str(tree_path))
        except Exception as e:
            pass  # Tree trimming failed, continue with node-based approach
    
    # Count events per clade pair
    pair_to_k = sub.groupby("clade_pair")["community_id"].nunique().to_dict()
    
    # Map each clade pair to MRCA node
    node_to_pairs = defaultdict(list)
    for (a,b), k in pair_to_k.items():
        mrca = find_mrca_node(a, b)
        if mrca is None:
            continue
        node_to_pairs[mrca].append(((a,b), k))
    
    # Apply parsimony rule at each node:
    # - If only one clade pair at this node: count its k
    # - If >1 distinct clade pairs: count 1 (but never below max k)
    comp_min = 0
    pernode = []
    for node, pairs in node_to_pairs.items():
        ks = [k for (_, k) in pairs]
        distinct_pairs = len(pairs)
        
        if distinct_pairs == 1:
            cnt = ks[0]
        else:
            # Multiple different clade pairs at same node = likely single ancestral transfer
            # But if one pair has multiple independent transfers (k>1), keep that
            cnt = max(1, max(ks))
        
        comp_min += cnt
        node_event_counts[node] += cnt
        pernode.append((node, distinct_pairs, max(ks), cnt))
    
    component_rows.append({
        "component_id": comp_id,
        "n_communities": sub["community_id"].nunique(),
        "n_clades": len(clades_in_component),
        "n_distinct_clade_pairs": len(pair_to_k),
        "min_events_component": comp_min
    })
    
    # Write per-component per-node breakdown (for debugging/analysis)
    if pernode:
        pd.DataFrame(pernode, columns=["node","n_pairs_at_node","max_k_samepair","events_counted"]).to_csv(
            out_events / f"{comp_id}_per_node.tsv", sep="\t", index=False
        )

# --- Write summaries
comp_df = pd.DataFrame(component_rows).sort_values("min_events_component", ascending=False)
comp_df.to_csv(out_summary / "minimal_events_by_component.tsv", sep="\t", index=False)

node_df = pd.DataFrame([
    {"node": n, "min_events_total": c} for n, c in node_event_counts.items()
]).sort_values("min_events_total", ascending=False)
node_df.to_csv(out_summary / "minimal_events_by_node.tsv", sep="\t", index=False)

total_min = int(comp_df["min_events_component"].sum()) if len(comp_df) else 0

summary_text = f"""Minimal HTT Events Summary
{'='*50}

Total components (single-linkage clusters): {len(comp_df)}
Total minimal HTT events: {total_min}

Events by node:
{node_df.to_string(index=False) if len(node_df) else 'No events'}

Component size distribution:
  Min communities per component: {comp_df['n_communities'].min() if len(comp_df) else 0}
  Max communities per component: {comp_df['n_communities'].max() if len(comp_df) else 0}
  Median: {comp_df['n_communities'].median() if len(comp_df) else 0}
"""

(out_summary / "minimal_events_summary.txt").write_text(summary_text)
(out_summary / "minimal_events_total.txt").write_text(f"Minimal HTT events (sum over components): {total_min}\n")

print("\nWrote:")
print(" -", out_summary / "minimal_events_by_component.tsv")
print(" -", out_summary / "minimal_events_by_node.tsv")
print(" -", out_summary / "minimal_events_summary.txt")
print(" -", out_summary / "minimal_events_total.txt")
print("")
print(f"TOTAL minimal HTT events: {total_min}")
PY

conda deactivate

echo ""
echo "### 08_minimal_htt_events finished: $(date)"
echo "Outputs:"
echo "  ${RESULTS}/network/community_to_component.tsv"
echo "  ${RESULTS}/network/component_sizes.tsv"
echo "  ${RESULTS}/trimmed_trees/  (if ete3 available)"
echo "  ${RESULTS}/events/  (per-component breakdowns)"
echo "  ${RESULTS}/summary/minimal_events_by_component.tsv"
echo "  ${RESULTS}/summary/minimal_events_by_node.tsv"
echo "  ${RESULTS}/summary/minimal_events_total.txt"