#!/bin/bash

### define path(s)
ROOT=/path/to/HTT
EARLGREY_RUNS=${ROOT}/results/04_results_te_annotation/earlgrey_runs
RESULTS=${ROOT}/results/04_results_te_annotation
LOGS=${ROOT}/logs/04_logs_te_annotation



### setup logging
JOBTAG="${SLURM_JOB_ID}_te_annotation_summary"
LOGFILE="${LOGS}/${JOBTAG}.log"
exec > >(tee "${LOGFILE}") 2>&1

echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"



### environment setup
export LANG=C
export LC_ALL=C
unset LANGUAGE

eval "$(conda shell.bash hook)"



### info
echo "Running on HYPERION - Date: $(date)"
echo "Results directory: ${RESULTS}"
echo "EarlGrey runs: ${EARLGREY_RUNS}"



### ===========================================================================
### Generate summary statistics
### ===========================================================================
echo ""
echo "### Generating summary statistics ###"
echo ""

conda activate HTT_python

python3 << EOF
from pathlib import Path

results_dir = Path("${RESULTS}")
earlgrey_runs = Path("${EARLGREY_RUNS}")
te_filtered_dir = results_dir / "te_filtered"
te_sequences_dir = results_dir / "te_sequences"
summary_dir = results_dir / "summary"
summary_dir.mkdir(exist_ok=True)

# collect statistics per genome
stats = []
total_tes = 0
total_raw = 0
total_filtered_out = 0

for gff_file in sorted(te_filtered_dir.glob("*_te_filtered.gff")):
    genome_name = gff_file.stem.replace("_te_filtered", "")
    
    # Find EarlGrey directory with glob pattern to handle version suffix
    earlgrey_dirs = list(earlgrey_runs.glob(f"{genome_name}_EarlGrey*"))
    raw_count = 0
    
    if earlgrey_dirs:
        earlgrey_dir = earlgrey_dirs[0]
        # Find the GFF file within the EarlGrey output
        raw_gff_list = list(earlgrey_dir.glob("*_summaryFiles/*.filteredRepeats.gff"))
        if raw_gff_list:
            raw_gff = raw_gff_list[0]
            with open(raw_gff) as f:
                for line in f:
                    if not line.startswith("#") and len(line.strip().split("\t")) >= 9:
                        raw_count += 1
    
    # count filtered TEs and classify
    te_classes = {}
    te_count = 0
    total_length = 0
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            te_count += 1
            length = int(fields[4]) - int(fields[3]) + 1
            total_length += length
            
            # TE class is in column 3 (GFF feature type)
            te_class = fields[2].split("/")[0]
            te_classes[te_class] = te_classes.get(te_class, 0) + 1
    
    filtered_out = raw_count - te_count
    pct_removed = (filtered_out / raw_count * 100) if raw_count > 0 else 0
    
    stats.append({
        'genome': genome_name,
        'raw_count': raw_count,
        'te_count': te_count,
        'filtered_out': filtered_out,
        'pct_removed': pct_removed,
        'total_bp': total_length,
        'classes': te_classes
    })
    total_tes += te_count
    total_raw += raw_count
    total_filtered_out += filtered_out
    
    print(f"{genome_name}: {raw_count} -> {te_count} TEs (removed {filtered_out}, {pct_removed:.1f}%), {total_length:,} bp")

# write summary table
summary_file = summary_dir / "te_annotation_summary.tsv"
with open(summary_file, "w") as out:
    out.write("genome\tte_count\ttotal_bp\n")
    for s in stats:
        out.write(f"{s['genome']}\t{s['te_count']}\t{s['total_bp']}\n")

print(f"\nTotal TEs across all genomes: {total_tes:,}")
print(f"Summary written to: {summary_file}")

# write filtering summary
filter_file = summary_dir / "te_filtering_summary.tsv"
with open(filter_file, "w") as out:
    out.write("genome\traw_count\tfiltered_count\tremoved_count\tpct_removed\n")
    for s in stats:
        out.write(f"{s['genome']}\t{s['raw_count']}\t{s['te_count']}\t{s['filtered_out']}\t{s['pct_removed']:.1f}\n")
    # totals row
    total_pct = (total_filtered_out / total_raw * 100) if total_raw > 0 else 0
    out.write(f"TOTAL\t{total_raw}\t{total_tes}\t{total_filtered_out}\t{total_pct:.1f}\n")

print(f"Filtering summary written to: {filter_file}")
print(f"\nOverall: {total_raw:,} -> {total_tes:,} TEs (removed {total_filtered_out:,}, {total_pct:.1f}%)")

# write class breakdown
class_file = summary_dir / "te_class_summary.tsv"
all_classes = set()
for s in stats:
    all_classes.update(s['classes'].keys())
all_classes = sorted(all_classes)

with open(class_file, "w") as out:
    out.write("genome\t" + "\t".join(all_classes) + "\n")
    for s in stats:
        row = [s['genome']] + [str(s['classes'].get(c, 0)) for c in all_classes]
        out.write("\t".join(row) + "\n")

print(f"Class breakdown written to: {class_file}")
EOF

conda deactivate

echo ""
echo "=========================================="
echo "Summary complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Summary statistics:   ${RESULTS}/summary/te_annotation_summary.tsv"
echo "  Filtering summary:    ${RESULTS}/summary/te_filtering_summary.tsv"
echo "  Class breakdown:      ${RESULTS}/summary/te_class_summary.tsv"
echo ""