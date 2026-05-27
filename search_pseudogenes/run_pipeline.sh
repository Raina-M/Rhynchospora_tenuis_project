#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Master runner for pseudogene detection pipeline
# =============================================================================
# USAGE:
#   bash run_pipeline.sh <miniprot.gff> <output_dir>
#
# OPTIONAL ENVIRONMENT VARIABLES:
#   FUNC_LENGTHS   path to functional gene span reference TSV
# =============================================================================

set -euo pipefail

GFF="${1:?Usage: bash run_pipeline.sh <miniprot.gff> <output_dir>}"
OUTDIR="${2:?Usage: bash run_pipeline.sh <miniprot.gff> <output_dir>}"
FUNC_LENGTHS="${FUNC_LENGTHS:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "========================================================="
echo " Pseudogene Detection Pipeline"
echo "========================================================="
echo " Input GFF:   $GFF"
echo " Output dir:  $OUTDIR"
echo ""

echo "[STEP 1/3] Parsing miniprot GFF..."
bash "$SCRIPT_DIR/01_parse_miniprot_gff.sh" "$GFF" "$OUTDIR"

echo ""
echo "[STEP 2/3] Classifying alignments..."
FUNC_ARG=""
if [[ -n "$FUNC_LENGTHS" ]]; then
    FUNC_ARG="--functional_gene_lengths $FUNC_LENGTHS"
fi

python3 "$SCRIPT_DIR/02_classify_pseudogenes.py" \
    --paf  "$OUTDIR/paf_table.tsv" \
    --mrna "$OUTDIR/mrna_table.tsv" \
    --cds  "$OUTDIR/cds_table.tsv" \
    --stop "$OUTDIR/stop_table.tsv" \
    $FUNC_ARG \
    --out  "$OUTDIR/pseudogene_candidates.tsv"

echo ""
echo "[STEP 3/3] Generating report and plots..."
Rscript "$SCRIPT_DIR/03_report_pseudogenes.R" \
    --input  "$OUTDIR/pseudogene_candidates.tsv" \
    --all    "$OUTDIR/pseudogene_candidates_all_alignments.tsv" \
    --outdir "$OUTDIR/report"

echo ""
echo "========================================================="
echo " Pipeline complete! Key outputs:"
echo "   $OUTDIR/pseudogene_candidates.tsv"
echo "   $OUTDIR/report/pseudogene_candidates.gff3"
echo "   $OUTDIR/report/pseudogene_diagnostic_plots.pdf"
echo "========================================================="
