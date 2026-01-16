#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
TREE_RESULTS=${ROOT}/results/01_results_phylo_tree
KS_RESULTS=${ROOT}/results/03_results_ks_divergence
CLUST_RESULTS=${ROOT}/results/07_results_htt_clustering
RESULTS=${ROOT}/results/08_results_minimal_htt_events
LOGS=${ROOT}/logs/08_logs_minimal_htt_events

### input files from previous pipeline steps
TREE_FILE=${TREE_RESULTS}/phylogeny_rooted.treefile
NODES_INFO=${KS_RESULTS}/node_analyses/nodes_info.json

### output from step 07
# Contains columns: community_id, qseqid, sseqid, qclade, sclade, ...
CLUSTERED_HTT_TSV=${CLUST_RESULTS}/summary/htt_pairs_with_community.tsv

### user-provided input
# IMPORTANT:
# Users MUST provide a clade assignment file mapping genomes to clades. The clade assigment file is generated using the local R scripts.
# This should be the same file used in previous steps.
CLADE_ASSIGN=${ROOT}/data/clade_assignment.tsv



### create results/log folder(s)
mkdir -p "${RESULTS}"/{network,clusters,events,trimmed_trees,summary,tmp}
mkdir -p "${LOGS}"



### setup logging
JOBTAG="${SLURM_JOB_ID}_minimal_htt_events"
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
echo "Tree file: ${TREE_FILE}"
echo "Nodes info: ${NODES_INFO}"
echo "Clade assignment: ${CLADE_ASSIGN}"
echo "Clustered HTT: ${CLUSTERED_HTT_TSV}"
echo "Results directory: ${RESULTS}"



### sanity checks
for f in "${TREE_FILE}" "${NODES_INFO}" "${CLADE_ASSIGN}" "${CLUSTERED_HTT_TSV}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        case "$f" in
            *01_results*) echo "Run 01_phylo_tree.sh first" ;;
            *03_results*) echo "Run 03_ks_divergence.sh first" ;;
            *07_results*) echo "Run 07_htt_clustering.sh first" ;;
            *clade_assignment*) echo "Please provide a clade assignment file" ;;
        esac
        exit 1
    fi
done

HTT_COUNT=$(tail -n +2 "${CLUSTERED_HTT_TSV}" | wc -l)
if [ "${HTT_COUNT}" -eq 0 ]; then
    echo "ERROR: No clustered HTT candidates found in ${CLUSTERED_HTT_TSV}"
    exit 1
fi
echo "Found ${HTT_COUNT} clustered HTT candidates"



### ===========================================================================
### STEP 1: Build network of communities and connected components
### ===========================================================================
echo ""
echo "### STEP 1: Building community network (single-linkage clustering) ###"
echo ""

conda deactivate
conda activate HTT_python

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"
echo "networkx version: $(python -c 'import networkx; print(networkx.__version__)')"

export CLUSTERED_HTT_TSV RESULTS

python3 << 'EOF'
import os
from pathlib import Path
import pandas as pd
import networkx as nx

clustered_tsv = Path(os.environ.get('CLUSTERED_HTT_TSV', ''))
results_dir = Path(os.environ.get('RESULTS', ''))
outdir = results_dir / "network"
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(clustered_tsv, sep="\t")
required = {"community_id","qseqid","sseqid","qclade","sclade"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"ERROR: {clustered_tsv} missing columns: {sorted(missing)}")

df["qseqid"] = df["qseqid"].astype(str)
df["sseqid"] = df["sseqid"].astype(str)
df["qclade"] = df["qclade"].astype(str)
df["sclade"] = df["sclade"].astype(str)
df["community_id"] = df["community_id"].astype(str)

# remove pairs where both strains belong to the same collapsed node
n_before = len(df)
collapsed_mask = df["qclade"].str.startswith("collapsed_node_")
same_clade_mask = df["qclade"] == df["sclade"]
df = df[~(collapsed_mask & same_clade_mask)].copy()
n_filtered = n_before - len(df)
print(f"Filtered {n_filtered} within-collapsed-node pairs")
print(f"Remaining HTT pairs: {len(df)}")

# build TE -> communities index
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

# single linkage clustering = connected components
components = list(nx.connected_components(G))
print(f"Connected components (single-linkage clusters): {len(components):,}")

# write mapping: community -> component
rows = []
for cid, comp in enumerate(components, start=1):
    for comm in comp:
        rows.append((comm, f"component_{cid:06d}", len(comp)))

map_df = pd.DataFrame(rows, columns=["community_id","component_id","component_n_communities"])
map_df.to_csv(outdir / "community_to_component.tsv", sep="\t", index=False)

# component summary
comp_sizes = map_df.groupby("component_id")["community_id"].count().reset_index(name="n_communities")
comp_sizes.to_csv(outdir / "component_sizes.tsv", sep="\t", index=False)

print(f"Wrote: {outdir / 'community_to_component.tsv'}")
print(f"Wrote: {outdir / 'component_sizes.tsv'}")
EOF

echo ""
echo "STEP 1 complete."

conda deactivate



### ===========================================================================
### STEP 2: Calculate minimal HTT events per component
### ===========================================================================
echo ""
echo "### STEP 2: Calculating minimal HTT events per component ###"
echo ""

conda activate HTT_python_ete3

echo "Running in environment: ${CONDA_DEFAULT_ENV}"
echo "Python version: $(python --version)"
echo "ete3 version: $(python -c 'import ete3; print(ete3.__version__)' 2>/dev/null || echo 'not available')"

export TREE_FILE NODES_INFO CLADE_ASSIGN CLUSTERED_HTT_TSV RESULTS

python3 << 'EOF'
import os
from pathlib import Path
import json
import pandas as pd
from collections import defaultdict

tree_file = Path(os.environ.get('TREE_FILE', ''))
nodes_info_file = Path(os.environ.get('NODES_INFO', ''))
clade_file = Path(os.environ.get('CLADE_ASSIGN', ''))
clustered_tsv = Path(os.environ.get('CLUSTERED_HTT_TSV', ''))
results_dir = Path(os.environ.get('RESULTS', ''))

map_file = results_dir / "network" / "community_to_component.tsv"
out_events = results_dir / "events"
out_trees = results_dir / "trimmed_trees"
out_summary = results_dir / "summary"
out_events.mkdir(parents=True, exist_ok=True)
out_trees.mkdir(parents=True, exist_ok=True)
out_summary.mkdir(parents=True, exist_ok=True)

# try to import ete3
try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False
    print("WARNING: ete3 not available. Tree trimming will be skipped.")

# load clade map: genome -> clade
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

# load species tree
if HAS_ETE3:
    tree_str = tree_file.read_text().strip()
    full_tree = Tree(tree_str, format=1)
    
    # simplify tip names to match clade_map keys
    for leaf in full_tree.iter_leaves():
        simple_name = leaf.name.split("_1_AssemblyScaffolds")[0]
        if simple_name.endswith("_1"):
            simple_name = simple_name[:-2]
        leaf.name = simple_name
    
    print(f"Loaded tree with {len(full_tree)} leaves")

# load nodes_info
nodes_info = json.loads(nodes_info_file.read_text())

# build node_to_clades sets and node_to_samples sets
node_to_clades = {}
node_to_samples = {}
node_size = {}
for node in nodes_info:
    clades = set()
    samples = set()
    for leaf in node["all_leaves"]:
        simple_leaf = leaf.split("_1_AssemblyScaffolds")[0]
        if simple_leaf.endswith("_1"):
            simple_leaf = simple_leaf[:-2]
        samples.add(simple_leaf)
        clades.add(clade_map.get(simple_leaf, simple_leaf))
    node_to_clades[node["node_name"]] = clades
    node_to_samples[node["node_name"]] = samples
    node_size[node["node_name"]] = len(clades)

# Identify collapsed clades and their member samples
collapsed_clades = {}
for sample, clade in clade_map.items():
    if clade.startswith("collapsed_node_"):
        if clade not in collapsed_clades:
            collapsed_clades[clade] = set()
        collapsed_clades[clade].add(sample)

print(f"Found {len(collapsed_clades)} collapsed clades:")
for cc, samples in collapsed_clades.items():
    print(f"  {cc}: {len(samples)} samples")

# Identify nodes that are internal to collapsed clades
# A node is internal to a collapsed clade if ALL its samples belong to that
# collapsed clade AND it has fewer samples than the collapsed clade
nodes_internal_to_collapsed = {}
for node_name, node_samples in node_to_samples.items():
    for cc_name, cc_samples in collapsed_clades.items():
        # Check if this node's samples are a proper subset of the collapsed clade
        if node_samples.issubset(cc_samples) and len(node_samples) < len(cc_samples):
            nodes_internal_to_collapsed[node_name] = cc_name
            break

print(f"Found {len(nodes_internal_to_collapsed)} nodes internal to collapsed clades")
for node, cc in nodes_internal_to_collapsed.items():
    print(f"  {node} is internal to {cc}")

def find_mrca_node(c1, c2):
    """Find MRCA node containing both clades, excluding internal collapsed nodes."""
    candidates = [
        n for n, cs in node_to_clades.items() 
        if c1 in cs and c2 in cs and n not in nodes_internal_to_collapsed
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda n: node_size[n])

# load HTT candidates with community assignments
df = pd.read_csv(clustered_tsv, sep="\t")
required = {"community_id","qseqid","sseqid","qclade","sclade"}
missing = required - set(df.columns)
if missing:
    raise SystemExit(f"ERROR: {clustered_tsv} missing columns: {sorted(missing)}")

df["qclade"] = df["qclade"].astype(str)
df["sclade"] = df["sclade"].astype(str)
df["community_id"] = df["community_id"].astype(str)

# remove pairs where both strains belong to the same collapsed node
n_before = len(df)
collapsed_mask = df["qclade"].str.startswith("collapsed_node_")
same_clade_mask = df["qclade"] == df["sclade"]
df = df[~(collapsed_mask & same_clade_mask)].copy()
n_filtered = n_before - len(df)
print(f"Filtered {n_filtered} within-collapsed-node pairs")
print(f"Remaining HTT pairs: {len(df)}")

# add component assignment
mp = pd.read_csv(map_file, sep="\t")
mp["community_id"] = mp["community_id"].astype(str)
df = df.merge(mp[["community_id","component_id"]], on="community_id", how="left")

if df["component_id"].isna().any():
    n_missing = df["component_id"].isna().sum()
    print(f"WARNING: {n_missing} rows without component assignment")
    df = df[df["component_id"].notna()].copy()

# deduplicate to unique community entries per clade pair
df_comm = df.drop_duplicates(subset=["community_id","qclade","sclade"]).copy()

def norm_pair(a,b):
    return tuple(sorted((a,b)))

df_comm["clade_pair"] = [norm_pair(a,b) for a,b in zip(df_comm["qclade"], df_comm["sclade"])]

# per component: trim tree and count minimal events
component_rows = []
node_event_counts = defaultdict(int)
total_components = df_comm["component_id"].nunique()

print(f"Processing {total_components} components...")

for comp_idx, (comp_id, sub) in enumerate(df_comm.groupby("component_id")):
    if (comp_idx + 1) % 100 == 0:
        print(f"  Processed {comp_idx + 1}/{total_components} components")
    
    # get all clades in this component
    clades_in_component = set(sub["qclade"]).union(set(sub["sclade"]))
    
    # get all genomes (samples) belonging to these clades
    samples_in_component = set()
    for sample, clade in clade_map.items():
        if clade in clades_in_component:
            samples_in_component.add(sample)
    
    # also add clades that are their own samples (uncollapsed)
    for clade in clades_in_component:
        if clade not in [c for c in clade_map.values()]:
            samples_in_component.add(clade)
    
    # trim tree to only these samples (if ete3 available)
    if HAS_ETE3 and len(samples_in_component) >= 2:
        try:
            leaves_to_keep = []
            for leaf in full_tree.iter_leaves():
                if leaf.name in samples_in_component:
                    leaves_to_keep.append(leaf.name)
            
            if len(leaves_to_keep) >= 2:
                trimmed = full_tree.copy()
                trimmed.prune(leaves_to_keep, preserve_branch_length=True)
                
                # save trimmed tree (only for manageable number of components)
                if total_components <= 1000:
                    tree_path = out_trees / f"{comp_id}.nwk"
                    trimmed.write(outfile=str(tree_path))
        except Exception:
            pass
    
    # count events per clade pair
    pair_to_k = sub.groupby("clade_pair")["community_id"].nunique().to_dict()
    
    # map each clade pair to MRCA node
    node_to_pairs = defaultdict(list)
    for (a,b), k in pair_to_k.items():
        mrca = find_mrca_node(a, b)
        if mrca is None:
            continue
        node_to_pairs[mrca].append(((a,b), k))
    
    # apply parsimony rule at each node
    comp_min = 0
    pernode = []
    for node, pairs in node_to_pairs.items():
        ks = [k for (_, k) in pairs]
        distinct_pairs = len(pairs)
        
        if distinct_pairs == 1:
            cnt = ks[0]
        else:
            # multiple different clade pairs at same node = likely single ancestral transfer
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
    
    # write per-component per-node breakdown
    if pernode:
        pd.DataFrame(pernode, columns=["node","n_pairs_at_node","max_k_samepair","events_counted"]).to_csv(
            out_events / f"{comp_id}_per_node.tsv", sep="\t", index=False
        )

# write summaries
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

print(f"\nWrote: {out_summary / 'minimal_events_by_component.tsv'}")
print(f"Wrote: {out_summary / 'minimal_events_by_node.tsv'}")
print(f"Wrote: {out_summary / 'minimal_events_summary.txt'}")
print(f"Wrote: {out_summary / 'minimal_events_total.txt'}")
print(f"\nTOTAL minimal HTT events: {total_min}")
EOF

echo ""
echo "STEP 2 complete."

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
echo "  Community to component:     ${RESULTS}/network/community_to_component.tsv"
echo "  Component sizes:            ${RESULTS}/network/component_sizes.tsv"
echo "  Trimmed trees:              ${RESULTS}/trimmed_trees/"
echo "  Per-component events:       ${RESULTS}/events/"
echo "  Events by component:        ${RESULTS}/summary/minimal_events_by_component.tsv"
echo "  Events by node:             ${RESULTS}/summary/minimal_events_by_node.tsv"
echo "  Summary:                    ${RESULTS}/summary/minimal_events_summary.txt"
echo "  Total events:               ${RESULTS}/summary/minimal_events_total.txt"
echo ""