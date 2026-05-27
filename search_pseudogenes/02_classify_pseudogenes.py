#!/usr/bin/env python3
"""
02_classify_pseudogenes.py
==============================================================================
PURPOSE:
    Classify each miniprot alignment as FUNCTIONAL, PROCESSED_PSEUDO,
    NONPROCESSED_PSEUDO, LOW_QUALITY_HIT, or AMBIGUOUS.

CORRECTIONS APPLIED:
  1. Relational Merge: Uses 'aln_idx' for explicit inner-joins instead of positional concat.
  2. Aligned Span Ratio: Bases 'span_ratio' strictly on the aligned query coordinates
     rather than total protein length to eliminate the truncation classification bug.
  3. Parent Gene Architecture: Added optional --query_exon_counts parameter to prevent 
     naturally single-exon orthologs from being falsely classified as PROCESSED_PSEUDO.
  4. Deactivated Unitary Rule: Masked genomes cannot reveal unitary losses safely. 
     Rule 4 has been deactivated to avoid blind-spot false positives.
  5. Logic Flow Repair: Rule 5 now tracks internal 'orf_disruptions', allowing pure 
     boundary truncations to gracefully land as clean partial-length FUNCTIONAL hits.
  6. String Defaults Guarded: Intercepts and parses "NA" or missing strings before 
     assigning boundary scores, removing the hidden 'da=0' (intact) corruption trap.
==============================================================================
"""

import argparse
import pandas as pd

# ==============================================================================
# THRESHOLDS
# ==============================================================================

MIN_QUERY_COVERAGE  = 0.50   # minimum fraction of protein that must be aligned
INTRON_RATIO_MAX    = 1.20   # genomic_span / (aligned_query_len*3) <= this → no introns
SINGLE_EXON_MAX_CDS = 2      # <=2 CDS blocks = "single exon" (processed pseudo signal)
INACTIVATING_MIN    = 1      # minimum total fs+st to call any pseudogene
MIN_IDENTITY        = 0.40   # below this → discard as non-homologous

# da/do thresholds for ORF boundary signals
DA_MISSING    = -1    # da=-1 means no upstream start codon found
DO_CLEAN      = 0     # do=0 means proper stop codon at C-terminus
DO_FAR_THRESH = 50    # do > this when q_coverage is high → likely missing stop codon

# ==============================================================================
# SAFE TYPE CONVERTERS
# ==============================================================================

def safe_int(value, default=0):
    """Convert to int with full fallback chain."""
    if value is None:
        return default
    try:
        if value != value:   # float NaN
            return default
    except TypeError:
        pass
    try:
        return int(float(str(value).strip()))
    except (ValueError, TypeError):
        pass
    parts = str(value).strip().split(":")
    if len(parts) == 3:
        try:
            return int(float(parts[2]))
        except (ValueError, TypeError):
            pass
    return default


def safe_float(value, default=0.0):
    """Convert to float with full fallback."""
    if value is None:
        return default
    try:
        if value != value:
            return default
    except TypeError:
        pass
    try:
        return float(str(value).strip())
    except (ValueError, TypeError):
        return default


# ==============================================================================
# CDS STATISTICS
# ==============================================================================

def compute_cds_stats(cds):
    """Per mRNA: count CDS blocks, total coding length, and splice site quality."""
    def summarise(g):
        n_blocks  = len(g)
        total_len = int((g["cds_end"] - g["cds_start"] + 1).sum())

        donors    = g["donor"].fillna("NA")
        acceptors = g["acceptor"].fillna("NA")
        has_noncanon = ((donors != "NA").any()) or ((acceptors != "NA").any())

        return pd.Series({
            "n_cds_blocks":      n_blocks,
            "total_cds_len":     total_len,
            "has_noncanon_splice": has_noncanon,
        })

    return cds.groupby("mrna_id").apply(summarise).reset_index()


# ==============================================================================
# STOP_CODON STATISTICS
# ==============================================================================

def compute_stop_stats(stop):
    """Count stop_codon GFF features per mRNA."""
    if stop is None or len(stop) == 0:
        return pd.DataFrame(columns=["mrna_id", "has_terminal_stop"])
    counts = (
        stop.groupby("mrna_id")
        .size()
        .rename("has_terminal_stop")
        .reset_index()
    )
    counts["has_terminal_stop"] = counts["has_terminal_stop"] > 0
    return counts


# ==============================================================================
# EVIDENCE COLLECTION
# ==============================================================================

def build_evidence(row, n_cds, has_noncanon_splice,
                   has_terminal_stop, is_only_hit, func_span, query_exon_count):
    """Collect all pseudogene evidence signals for one alignment."""
    ev = []

    # ── Premature frameshifts and stop codons ─────────────────────────────────
    fs     = safe_int(row.get("frameshifts_paf", 0))
    st     = safe_int(row.get("stop_codons_paf", 0))
    fs_gff = safe_int(row.get("frameshift_gff",  0))
    st_gff = safe_int(row.get("stopcodon_gff",   0))

    fs_total = max(fs, fs_gff)
    st_total = max(st, st_gff)

    if fs_total > 0:
        ev.append(f"FRAMESHIFT:fs_paf={fs};fs_gff={fs_gff}")
    if st_total > 0:
        ev.append(f"PREMATURE_STOP:st_paf={st};st_gff={st_gff}")

    # ── da: distance to nearest upstream START codon ──────────────────────────
    raw_da = str(row.get("da", "NA")).strip()
    if raw_da in ["NA", "", "nan"]:
        da = 999
    else:
        da = safe_int(raw_da, default=999)

    if raw_da == "-1" or da == -1:
        ev.append("MISSING_START:da=-1(no_upstream_ATG)")
    elif da == 0:
        ev.append("INTACT_START:da=0(ATG_at_alignment_start)")
    else:
        ev.append(f"PARTIAL_START:da={da}bp_upstream")

    # ── do: distance to nearest downstream STOP codon ────────────────────────
    raw_do = str(row.get("do", "NA")).strip()
    if raw_do in ["NA", "", "nan"]:
        do = 999
    else:
        do = safe_int(raw_do, default=999)

    q_cov = safe_float(row.get("q_coverage", 0.0))
    if do == 0:
        ev.append("INTACT_STOP:do=0(stop_codon_at_alignment_end)")
    elif do > DO_FAR_THRESH and q_cov >= 0.90:
        ev.append(f"MISSING_STOP:do={do}bp;q_cov={q_cov:.2f}(stop_codon_likely_lost)")
    else:
        ev.append(f"DISTANT_STOP:do={do}bp")

    # ── Terminal stop_codon GFF feature ───────────────────────────────────────
    if has_terminal_stop:
        ev.append("TERMINAL_STOP_FEATURE:proper_C-terminal_termination")

    # ── Exon / intron structure ───────────────────────────────────────────────
    if n_cds <= SINGLE_EXON_MAX_CDS:
        ev.append(f"SINGLE_EXON:n_cds={n_cds}")
    else:
        ev.append(f"MULTI_EXON:n_cds={n_cds}")

    if query_exon_count is not None:
        ev.append(f"QUERY_PARENT_EXONS:{query_exon_count}")

    # ── Splice site canonicality ──────────────────────────────────────────────
    if has_noncanon_splice:
        ev.append("NONCANON_SPLICE:non-GT-AG_intron_detected")
    elif n_cds > SINGLE_EXON_MAX_CDS:
        ev.append("CANONICAL_SPLICE:GT-AG_inferred(no_noncanon_attrs)")

    # ── Aligned Genomic span ratio (Truncation-Aware Fix) ─────────────────────
    q_start = safe_int(row.get("q_start", 0))
    q_end = safe_int(row.get("q_end", 0))
    aligned_query_len = max(1, q_end - q_start)
    expected_mrna_len = aligned_query_len * 3
    genomic_span      = safe_int(row.get("genomic_span", 0))
    span_ratio        = genomic_span / expected_mrna_len

    ev.append(f"SPAN_RATIO:{span_ratio:.2f}(genomic_bp/aligned_query_mRNA_bp)")
    if span_ratio <= INTRON_RATIO_MAX:
        ev.append("COMPACT_SPAN:no_intron_space")

    if func_span is not None and func_span > 0:
        ev.append(f"SPAN_VS_FUNCTIONAL:{genomic_span/func_span:.2f}")

    # ── Alignment quality ─────────────────────────────────────────────────────
    identity = safe_float(row.get("identity", 0.0))
    ev.append(f"QUERY_COVERAGE:{q_cov:.3f}")
    ev.append(f"IDENTITY:{identity:.3f}")
    if q_cov < 0.80:
        ev.append(f"TRUNCATED_ALIGNMENT:q_cov={q_cov:.3f}")

    # ── Rank and uniqueness ───────────────────────────────────────────────────
    rank = safe_int(row.get("rank", 1))
    ev.append(f"RANK:{rank}")
    if is_only_hit:
        ev.append("ONLY_GENOMIC_HIT:one_locus_for_this_protein")

    # ── mapq ──────────────────────────────────────────────────────────────────
    mapq = safe_int(row.get("mapq", 0))
    if mapq == 255:
        ev.append("MAPQ_MISSING:255")

    return ev, fs_total, st_total


# ==============================================================================
# CLASSIFICATION
# ==============================================================================

def classify(row, n_cds, evidence, fs_total, st_total,
             has_noncanon_splice, has_terminal_stop, is_only_hit, query_exon_count):
    """Assign classification and confidence score."""
    identity     = safe_float(row.get("identity",   0.0))
    q_cov        = safe_float(row.get("q_coverage", 0.0))
    genomic_span = safe_int(row.get("genomic_span", 0))
    
    q_start = safe_int(row.get("q_start", 0))
    q_end = safe_int(row.get("q_end", 0))
    aligned_query_len = max(1, q_end - q_start)
    span_ratio   = genomic_span / (aligned_query_len * 3)

    # da/do signals
    raw_da = str(row.get("da", "NA")).strip()
    if raw_da in ["NA", "", "nan"]:
        da = 999
    else:
        da = safe_int(raw_da, default=999)

    raw_do = str(row.get("do", "NA")).strip()
    if raw_do in ["NA", "", "nan"]:
        do = 999
    else:
        do = safe_int(raw_do, default=999)

    missing_start = (raw_da == "-1" or da == -1)
    missing_stop  = (do > DO_FAR_THRESH and q_cov >= 0.90)

    # Count ORF disruption signals
    orf_disruptions = fs_total + st_total
    boundary_disruptions = sum([missing_start, missing_stop])
    total_disruptions = orf_disruptions + boundary_disruptions

    is_compact    = span_ratio <= INTRON_RATIO_MAX
    is_single_exon= n_cds <= SINGLE_EXON_MAX_CDS

    # Guard against naturally single-exon parent genes
    is_parent_multi_exon = True
    if query_exon_count is not None and query_exon_count <= SINGLE_EXON_MAX_CDS:
        is_parent_multi_exon = False

    # Rule 0: discard low-quality hits
    if identity < MIN_IDENTITY or q_cov < MIN_QUERY_COVERAGE:
        return "LOW_QUALITY_HIT", 0.0

    # Rule 1: no disruptions at all → functional
    if total_disruptions == 0:
        return "FUNCTIONAL", min(1.0, identity * q_cov)

    # Rule 2: processed pseudogene
    proc_signals = sum([is_compact, is_single_exon])
    if orf_disruptions >= INACTIVATING_MIN and proc_signals == 2 and is_parent_multi_exon:
        conf = min(1.0, (proc_signals / 2.0) * min(1.0, total_disruptions * 0.3 + 0.5))
        return "PROCESSED_PSEUDO", conf
    if total_disruptions >= INACTIVATING_MIN and proc_signals == 2 and is_parent_multi_exon:
        conf = min(1.0, 0.55 + total_disruptions * 0.1)
        return "PROCESSED_PSEUDO", conf

    # Rule 3: non-processed pseudogene (includes single-exon parents undergoing standard duplication)
    nonproc_signals = sum([not is_single_exon, not is_compact])
    if orf_disruptions >= INACTIVATING_MIN and (nonproc_signals >= 1 or not is_parent_multi_exon):
        conf = min(1.0, (max(1, nonproc_signals) / 2.0 + 0.33) * min(1.0, total_disruptions * 0.3 + 0.4))
        return "NONPROCESSED_PSEUDO", conf

    # Rule 4: unitary pseudogene -> DEACTIVATED due to genome masking
    # True unitary status can only be reliably judged on clean unmasked genomes.
    # if is_only_hit and total_disruptions >= INACTIVATING_MIN:
    #     return "UNITARY_PSEUDO", min(1.0, 0.5 + total_disruptions * 0.12)

    # Rule 5: ambiguous structure (has internal ORF disruptions but ambiguous structures)
    if orf_disruptions >= INACTIVATING_MIN:
        return "AMBIGUOUS", 0.3

    # Clean partial/edge truncations without active internal disruptions drop out as FUNCTIONAL
    return "FUNCTIONAL", min(1.0, identity * q_cov)


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Classify miniprot alignments as pseudogene types.")
    parser.add_argument("--paf",  required=True)
    parser.add_argument("--mrna", required=True)
    parser.add_argument("--cds",  required=True)
    parser.add_argument("--stop", required=True)
    parser.add_argument("--functional_gene_lengths", default=None)
    parser.add_argument("--query_exon_counts", default=None, 
                        help="Optional tab-separated file mapping protein_id to natural exon count.")
    parser.add_argument("--out",  required=True)
    args = parser.parse_args()

    print("[02] Loading tables...")
    paf  = pd.read_csv(args.paf,  sep="\t")
    mrna = pd.read_csv(args.mrna, sep="\t")
    cds  = pd.read_csv(args.cds,  sep="\t")
    stop = pd.read_csv(args.stop, sep="\t")
    print(f"     PAF:{len(paf)} | mRNA:{len(mrna)} | CDS:{len(cds)} | stop_codon:{len(stop)}")

    func_spans = {}
    if args.functional_gene_lengths:
        fl = pd.read_csv(args.functional_gene_lengths, sep="\t",
                         header=None, names=["protein_id", "func_span"])
        func_spans = dict(zip(fl["protein_id"].astype(str), fl["func_span"].astype(float)))
        print(f"     Functional spans: {len(func_spans)} proteins")

    query_exons = {}
    if args.query_exon_counts:
        qe = pd.read_csv(args.query_exon_counts, sep="\t", header=None, names=["protein_id", "exon_count"])
        query_exons = dict(zip(qe["protein_id"].astype(str), qe["exon_count"].astype(int)))
        print(f"     Query exon counts: {len(query_exons)} proteins")

    print("[02] Computing CDS and stop_codon stats...")
    cds_stats  = compute_cds_stats(cds)
    stop_stats = compute_stop_stats(stop)

    # Relational merge via aln_idx to preserve order and structure safely
    merged = pd.merge(
        paf,
        mrna[["aln_idx", "mrna_id", "rank", "identity", "positive",
              "frameshift_gff", "stopcodon_gff"]],
        on="aln_idx",
        how="inner"
    )

    merged = merged.merge(cds_stats, on="mrna_id", how="left")
    merged["n_cds_blocks"]       = merged["n_cds_blocks"].fillna(0).astype(int)
    merged["total_cds_len"]      = merged["total_cds_len"].fillna(0).astype(int)
    merged["has_noncanon_splice"]= merged["has_noncanon_splice"].fillna(False)

    merged = merged.merge(stop_stats, on="mrna_id", how="left")
    merged["has_terminal_stop"]  = merged["has_terminal_stop"].fillna(False)

    hit_counts = (
        merged.groupby("protein_id")["mrna_id"]
        .count().rename("total_hits").reset_index()
    )
    merged = merged.merge(hit_counts, on="protein_id")
    merged["is_only_hit"] = merged["total_hits"] == 1

    print("[02] Classifying alignments...")
    results = []

    for _, row in merged.iterrows():
        n_cds         = int(row["n_cds_blocks"])
        has_noncanon  = bool(row["has_noncanon_splice"])
        has_term_stop = bool(row["has_terminal_stop"])
        is_only       = bool(row["is_only_hit"])
        func_span     = func_spans.get(str(row["protein_id"]), None)
        query_exon_count = query_exons.get(str(row["protein_id"]), None)

        evidence, fs_total, st_total = build_evidence(
            row, n_cds, has_noncanon, has_term_stop, is_only, func_span, query_exon_count)
        label, confidence = classify(
            row, n_cds, evidence, fs_total, st_total,
            has_noncanon, has_term_stop, is_only, query_exon_count)

        prot_len     = max(1, safe_int(row.get("protein_len")))
        genomic_span = safe_int(row.get("genomic_span"))
        
        q_start = safe_int(row.get("q_start", 0))
        q_end = safe_int(row.get("q_end", 0))
        aligned_query_len = max(1, q_end - q_start)

        da_val       = str(row.get("da", "NA")).strip()
        do_val       = str(row.get("do", "NA")).strip()

        results.append({
            "mrna_id":                  row["mrna_id"],
            "protein_id":               row["protein_id"],
            "contig":                   row.get("contig", "NA"),
            "g_start":                  row.get("g_start", "NA"),
            "g_end":                    row.get("g_end",   "NA"),
            "strand":                   row.get("strand",  "NA"),
            "protein_len_aa":           prot_len,
            "q_coverage":               round(safe_float(row.get("q_coverage")), 3),
            "identity":                 round(safe_float(row.get("identity")),   3),
            "mapq":                     safe_int(row.get("mapq")),
            "rank":                     safe_int(row.get("rank"), default=1),
            "genomic_span_bp":          genomic_span,
            "expected_mrna_bp":         aligned_query_len * 3,
            "span_ratio":               round(genomic_span / (aligned_query_len * 3), 3),
            "frameshifts":              fs_total,
            "premature_stop_codons":    st_total,
            "da":                       da_val,
            "do":                       do_val,
            "missing_start_codon":      (da_val == "-1"),
            "clean_stop_codon":         (do_val == "0"),
            "n_cds_blocks":             n_cds,
            "has_noncanon_splice":      has_noncanon,
            "has_terminal_stop_feat":   has_term_stop,
            "total_hits_this_protein":  int(row["total_hits"]),
            "is_only_hit":              is_only,
            "classification":           label,
            "confidence":               round(confidence, 3),
            "evidence":                 " | ".join(evidence),
        })

    out_df = pd.DataFrame(results)

    print("\n[02] Classification summary:")
    for label, count in out_df["classification"].value_counts().items():
        print(f"     {label:<25s}  {count:>6d}")

    all_out = args.out.replace(".tsv", "_all_alignments.tsv")
    out_df.to_csv(all_out, sep="\t", index=False)
    print(f"\n[02] All alignments        → {all_out}")

    pseudo_labels = {"PROCESSED_PSEUDO", "NONPROCESSED_PSEUDO",
                     "UNITARY_PSEUDO", "AMBIGUOUS"}
    pseudo_df = (
        out_df[out_df["classification"].isin(pseudo_labels)]
        .copy()
        .sort_values(["classification", "confidence"], ascending=[True, False])
    )
    pseudo_df.to_csv(args.out, sep="\t", index=False)
    print(f"[02] Pseudogene candidates → {args.out}  ({len(pseudo_df)} rows)")


if __name__ == "__main__":
    main()