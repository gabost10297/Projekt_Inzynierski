#!/bin/bash
set -euo pipefail

# Config
BLAST_DB="/data/database/unite_blast_db"
THREADS="${BLAST_THREADS:-10}"
MAX_TARGET_SEQS=5
PERC_IDENTITY=90
EVALUE_CUTOFF=1e-20
# Confidence rating thresholds (percent identity / query coverage)
HIGH_PIDENT=99
HIGH_QCOV=92
LOW_PIDENT=97
LOW_QCOV=85

mkdir -p /data/blast_results

echo "STARTING LOCAL BLAST SEARCH"
echo "====================================================="

for sample_dir in /data/consensus_results/*; do
    [[ -d "$sample_dir" ]] || continue

    SAMPLE=$(basename "$sample_dir")
    echo ">>> BLASTing sample: ${SAMPLE} <<<"

    OUT_ALL="/data/blast_results/${SAMPLE}_blast_summary.tsv"
    OUT_STRICT="/data/blast_results/${SAMPLE}_blast_summary_strict.tsv"

    HEADER=$'Cluster_Name\tQuery_Length\tLength_Tier\tReference_ID\tPercent_Identity\tAlignment_Length\tQuery_Coverage(%)\tE-value\tBitscore\tTop2_Pident_Gap\tSpecies_Name\tConfidence\tAssigned_Level'
    echo -e "$HEADER" > "$OUT_ALL"

    for fasta in "$sample_dir"/*.fasta; do
        [[ -f "$fasta" ]] || continue

        CLUSTER=$(basename "$fasta" .fasta)

        QUERY_LEN=$(grep -v '^>' "$fasta" | tr -d '\n' | wc -c)
        if (( QUERY_LEN < 350 )); then
            LENGTH_TIER="D"
        elif (( QUERY_LEN < 400 )); then
            LENGTH_TIER="C"
        elif (( QUERY_LEN < 500 )); then
            LENGTH_TIER="B"
        elif (( QUERY_LEN <= 900 )); then
            LENGTH_TIER="A"
        else
            LENGTH_TIER="LONG"
        fi

        TMP_BLAST=$(mktemp)
        if ! blastn -query "$fasta" -db "$BLAST_DB" \
            -task blastn \
            -perc_identity "$PERC_IDENTITY" \
            -evalue "$EVALUE_CUTOFF" \
            -max_target_seqs "$MAX_TARGET_SEQS" \
            -num_threads "$THREADS" \
            -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle" \
            -out "$TMP_BLAST" 2>/dev/null; then
            :
        fi

        if [[ ! -s "$TMP_BLAST" ]]; then
            echo -e "${CLUSTER}\t${QUERY_LEN}\t${LENGTH_TIER}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tfail\tunassigned" >> "$OUT_ALL"
            rm -f "$TMP_BLAST"
            continue
        fi

        awk -v cluster="$CLUSTER" -v qlen="$QUERY_LEN" -v tier="$LENGTH_TIER" \
            -v high_pid="$HIGH_PIDENT" -v high_qcov="$HIGH_QCOV" \
            -v low_pid="$LOW_PIDENT" -v low_qcov="$LOW_QCOV" '
        BEGIN { FS = OFS = "\t" }
        NR == 1 {
            ref = $2
            p1 = $3 + 0
            aln = $4 + 0
            qcov = $5 + 0
            eval = $6
            bit = $7 + 0
            title = $8
            for (i = 9; i <= NF; i++) title = title " " $i
        }
        NR == 2 { p2 = $3 + 0 }
        END {
            gap = (NR >= 2 ? p1 - p2 : 100)
            if (gap < 0) gap = 0

            conf = "fail"
            level = "unassigned"

            if (tier == "LONG") {
                conf = "review_long"
                level = "manual_review"
            } else if (tier == "A" && p1 >= high_pid && qcov >= high_qcov && aln >= 400 && eval <= 1e-20 && gap >= 0.5) {
                conf = "high"
                level = "species"
            } else if (tier == "A" && p1 >= low_pid && qcov >= low_qcov && aln >= 400 && eval <= 1e-20 && gap >= 0.5) {
                conf = "medium"
                level = "species"
            } else if (tier == "B" && p1 >= high_pid && qcov >= high_qcov && aln >= 350 && gap >= 0.5) {
                conf = "high"
                level = "species"
            } else if (tier == "B" && p1 >= low_pid && qcov >= low_qcov && aln >= 350 && gap >= 0.5) {
                conf = "medium"
                level = "species"
            } else if (tier == "C" && p1 >= low_pid && qcov >= low_qcov && gap >= 1.0) {
                conf = "low_species"
                level = "species"
            } else if (p1 >= 90 && qcov >= 75 && aln >= 200) {
                conf = "genus_only"
                level = "genus"
            }

            if (tier == "D" && level == "species") {
                conf = "fail"
                level = "genus"
            }
            if (gap < 0.5 && level == "species") {
                conf = "ambiguous"
                level = "genus"
            }

            print cluster, qlen, tier, ref, p1, aln, qcov, eval, bit, gap, title, conf, level
        }' "$TMP_BLAST" >> "$OUT_ALL"

        rm -f "$TMP_BLAST"
    done

    awk -F'\t' -v OFS='\t' 'NR==1 { print; next } $(NF-1)=="high" || $(NF-1)=="medium" { print }' \
        "$OUT_ALL" > "$OUT_STRICT"

    echo "Done: ${SAMPLE}"
    echo "  All:    ${OUT_ALL}"
    echo "  Strict: ${OUT_STRICT}"
done

echo "====================================================="
echo "ALL SAMPLES IDENTIFIED!"