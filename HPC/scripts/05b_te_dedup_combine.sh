#!/bin/bash
#SBATCH --job-name=te_dedup_combine
#SBATCH --qos hprio
#SBATCH --account node
#SBATCH --partition node
#SBATCH --mail-user=felix.zimmermann@wsl.ch
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=/storage/zimmermf/HTT/logs/05_logs_te_classification/dedup_combine_%j.out
#SBATCH --error=/storage/zimmermf/HTT/logs/05_logs_te_classification/dedup_combine_%j.err

set -euo pipefail

ROOT=/storage/zimmermf/HTT
RESULTS=${ROOT}/results/05_results_te_classification
FILTERED_DIR=${RESULTS}/filtered_tes
DEDUP_DIR=${RESULTS}/deduplicated
SUMMARY_DIR=${RESULTS}/summary

mkdir -p "${DEDUP_DIR}" "${SUMMARY_DIR}"

eval "$(conda shell.bash hook)"

echo "Running on HYPERION - Date: $(date)"
echo "Results directory: ${RESULTS}"

### ===========================================================================
### STEP 1: Deduplicate TEs per genome (exact identical sequences)
### STEP 2: Combine all deduplicated TEs
### ===========================================================================
echo ""
echo "### Deduplicating and combining TEs ###"
echo ""

conda activate HTT_python

python3 - "${FILTERED_DIR}" "${DEDUP_DIR}" "${RESULTS}" << 'EOF'
import sys
from pathlib import Path
from hashlib import md5

filtered_dir = Path(sys.argv[1])
dedup_dir = Path(sys.argv[2])
results_dir = Path(sys.argv[3])

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
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

# Process each strain
stats = []
all_tes = []

for fasta in sorted(filtered_dir.glob("*_filtered_tes.fasta")):
    strain = fasta.stem.replace("_filtered_tes", "")
    
    # Read all sequences
    sequences = list(read_fasta(fasta))
    before_count = len(sequences)
    
    # Deduplicate by exact sequence hash
    seen_hashes = {}
    dedup_records = []
    
    for header, seq in sequences:
        seq_hash = md5(seq.encode()).hexdigest()
        if seq_hash not in seen_hashes:
            seen_hashes[seq_hash] = header
            dedup_records.append((header, seq))
    
    after_count = len(dedup_records)
    removed = before_count - after_count
    
    # Write deduplicated fasta
    outfile = dedup_dir / f"{strain}_dedup.fasta"
    write_fasta(outfile, dedup_records)
    
    # Add to combined list with strain prefix
    for header, seq in dedup_records:
        new_id = f"{strain}|{header}"
        all_tes.append((new_id, strain, header, seq))
    
    stats.append({
        "strain": strain,
        "before": before_count,
        "after": after_count,
        "removed": removed
    })
    
    print(f"{strain}: {before_count} -> {after_count} TEs (removed {removed} identical)")

# Write combined fasta
combined_file = results_dir / "all_tes_combined.fasta"
with open(combined_file, "w") as f:
    for new_id, strain, orig_id, seq in all_tes:
        f.write(f">{new_id}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

# Write info table
info_file = results_dir / "all_tes_info.tsv"
with open(info_file, "w") as f:
    f.write("te_id\tstrain\toriginal_id\n")
    for new_id, strain, orig_id, seq in all_tes:
        f.write(f"{new_id}\t{strain}\t{orig_id}\n")

# Write dedup summary
summary_file = results_dir / "summary" / "dedup_summary.tsv"
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

print(f"\n{'='*50}")
print(f"Total: {total_before:,} -> {total_after:,} TEs")
print(f"Removed {total_removed:,} identical sequences ({total_removed/total_before*100:.1f}%)")
print(f"{'='*50}")
print(f"\nOutput files:")
print(f"  Combined FASTA: {combined_file}")
print(f"  TE info table:  {info_file}")
print(f"  Dedup summary:  {summary_file}")
EOF

conda deactivate

echo ""
echo "Done!"