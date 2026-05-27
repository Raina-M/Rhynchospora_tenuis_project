#!/usr/bin/env bash
# =============================================================================
# 01_parse_miniprot_gff.sh
# =============================================================================
# PURPOSE:
#   Parse a miniprot GFF3 file into four clean TSV tables:
#     paf_table.tsv   — one row per aligned protein (from ##PAF lines)
#     mrna_table.tsv  — one row per mRNA GFF feature
#     cds_table.tsv   — one row per CDS GFF feature
#     stop_table.tsv  — one row per stop_codon GFF feature
#
# CORRECTIONS APPLIED:
#   - Added an incremental alignment index (aln_idx) to both PAF and mRNA 
#     tables to enable a robust relational merge instead of a risky positional merge.
# =============================================================================

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <miniprot.gff> <output_dir>"
    exit 1
fi

GFF="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

echo "[01] Parsing miniprot GFF: $GFF"
echo "[01] Output directory:     $OUTDIR"

# =============================================================================
# PART 1 — ##PAF lines
# =============================================================================

echo "[01] Extracting PAF lines..."

echo -e "aln_idx\tprotein_id\tprotein_len\tq_start\tq_end\tq_coverage\tstrand\tcontig\tg_start\tg_end\tgenomic_span\tn_match_nt\taln_span_excl_intron\tmapq\tframeshifts_paf\tstop_codons_paf\tda\tdo\traw_tags" \
    > "$OUTDIR/paf_table.tsv"

grep "^##PAF" "$GFF" | awk '
BEGIN { OFS="\t"; aln_idx=0 }
{
    # Skip unaligned proteins (strand field is "*")
    if ($6 == "*") next

    aln_idx++
    protein_id          = $2
    protein_len         = $3
    q_start             = $4    # 0-based start on protein
    q_end               = $5    # 0-based end on protein
    strand              = $6
    contig              = $7
    # $8 = contig_len (skip)
    g_start             = $9    # 0-based genomic start
    g_end               = $10   # 0-based genomic end
    n_match_nt          = $11   # matching nt (introns excluded)
    aln_span_excl_intron= $12   # alignment block size (introns excluded)
    mapq                = $13   # 0-255; 255 = missing

    q_coverage   = (protein_len > 0) ? (q_end - q_start) / protein_len : 0
    # Genomic span INCLUDES introns — use this for intron detection
    genomic_span = g_end - g_start

    # Scan optional tags: KEY:TYPE:VALUE (semicolons in raw_tags for TSV safety)
    fs       = 0
    st       = 0
    da       = "NA"
    do_val   = "NA"
    raw_tags = ""

    for (i = 14; i <= NF; i++) {
        tag = $i
        raw_tags = (raw_tags == "") ? tag : raw_tags ";" tag
        n = split(tag, a, ":")
        if (n < 3) continue
        key = a[1]; val = a[3]
        if (key == "fs") fs     = int(val)
        if (key == "st") st     = int(val)
        if (key == "da") da     = val   # keep as string: -1 / 0 / positive int
        if (key == "do") do_val = val   # keep as string: -1 / 0 / positive int
    }

    print aln_idx, protein_id, protein_len, q_start, q_end, q_coverage, \
          strand, contig, g_start, g_end, genomic_span, \
          n_match_nt, aln_span_excl_intron, mapq, \
          fs, st, da, do_val, raw_tags
}' >> "$OUTDIR/paf_table.tsv"

echo "[01] PAF table: $(( $(wc -l < "$OUTDIR/paf_table.tsv") - 1 )) aligned proteins"

# =============================================================================
# PART 2 — mRNA feature lines
# =============================================================================

echo "[01] Extracting mRNA features..."

echo -e "aln_idx\tmrna_id\tprotein_id\tcontig\tg_start\tg_end\tstrand\tscore\trank\tidentity\tpositive\tframeshift_gff\tstopcodon_gff\tq_start_target\tq_end_target" \
    > "$OUTDIR/mrna_table.tsv"

awk '
BEGIN { OFS="\t"; aln_idx=0 }
$3 == "mRNA" {
    aln_idx++
    contig  = $1; g_start = $4; g_end = $5
    score   = $6; strand  = $7; attrs = $9

    mrna_id   = extract(attrs, "ID")
    rank      = extract(attrs, "Rank")
    identity  = extract(attrs, "Identity")
    positive  = extract(attrs, "Positive")
    fs_gff    = extract(attrs, "Frameshift")   # NA if absent
    st_gff    = extract(attrs, "StopCodon")    # NA if absent

    # Target= value is: "protein_id q_start q_end" (space-separated)
    target_raw = extract(attrs, "Target")
    n = split(target_raw, t, " ")
    protein_id = (n >= 1) ? t[1] : "NA"
    q_start_t  = (n >= 2) ? t[2] : "NA"
    q_end_t    = (n >= 3) ? t[3] : "NA"

    print aln_idx, mrna_id, protein_id, contig, g_start, g_end, strand, score, \
          rank, identity, positive, fs_gff, st_gff, q_start_t, q_end_t
}

function extract(attrs, key,    n, parts, i) {
    n = split(attrs, parts, ";")
    for (i = 1; i <= n; i++) {
        if (parts[i] ~ "^" key "=") {
            sub("^" key "=", "", parts[i])
            return parts[i]
        }
    }
    return "NA"
}
' "$GFF" >> "$OUTDIR/mrna_table.tsv"

echo "[01] mRNA table: $(( $(wc -l < "$OUTDIR/mrna_table.tsv") - 1 )) mRNA records"

# =============================================================================
# PART 3 — CDS feature lines
# =============================================================================

echo "[01] Extracting CDS features..."

echo -e "mrna_id\tcontig\tcds_start\tcds_end\tstrand\tphase\tcds_identity\tframeshift_cds\tstopcodon_cds\tdonor\tacceptor" \
    > "$OUTDIR/cds_table.tsv"

awk '
$3 == "CDS" {
    contig    = $1; cds_start = $4; cds_end = $5
    strand    = $7; phase     = $8; attrs   = $9

    parent       = extract(attrs, "Parent")
    cds_identity = extract(attrs, "Identity")
    fs_cds       = extract(attrs, "Frameshift")
    st_cds       = extract(attrs, "StopCodon")
    donor        = extract(attrs, "Donor")
    acceptor     = extract(attrs, "Acceptor")

    print parent, contig, cds_start, cds_end, strand, phase, \
          cds_identity, fs_cds, st_cds, donor, acceptor
}

function extract(attrs, key,    n, parts, i) {
    n = split(attrs, parts, ";")
    for (i = 1; i <= n; i++) {
        if (parts[i] ~ "^" key "=") {
            sub("^" key "=", "", parts[i])
            return parts[i]
        }
    }
    return "NA"
}
' OFS="\t" "$GFF" >> "$OUTDIR/cds_table.tsv"

echo "[01] CDS table: $(( $(wc -l < "$OUTDIR/cds_table.tsv") - 1 )) CDS records"

# =============================================================================
# PART 4 — stop_codon feature lines
# =============================================================================

echo "[01] Extracting stop_codon features..."

echo -e "mrna_id\tcontig\tstart\tend\tstrand" > "$OUTDIR/stop_table.tsv"

awk '
$3 == "stop_codon" {
    contig = $1; start = $4; end = $5; strand = $7; attrs = $9
    parent = extract(attrs, "Parent")
    print parent, contig, start, end, strand
}
function extract(attrs, key,    n, parts, i) {
    n = split(attrs, parts, ";")
    for (i = 1; i <= n; i++) {
        if (parts[i] ~ "^" key "=") {
            sub("^" key "=", "", parts[i])
            return parts[i]
        }
    }
    return "NA"
}
' OFS="\t" "$GFF" >> "$OUTDIR/stop_table.tsv"

echo "[01] stop_codon table: $(( $(wc -l < "$OUTDIR/stop_table.tsv") - 1 )) records"
echo ""
echo "[01] Done."
ls -lh "$OUTDIR"