#!/usr/bin/env Rscript
# ==============================================================================
# 03_report_pseudogenes.R
# ==============================================================================
# PURPOSE:
#   Read pseudogene candidate TSVs and produce:
#     1. Structured text report (stdout)
#     2. Diagnostic plots PDF (8 panels including premature stop codon stats)
#     3. Clean GFF3 output file
#
# CORRECTIONS from reading the miniprot manual:
#   - stop_codon GFF feature = proper C-terminal termination (GOOD sign).
#     Renamed to has_terminal_stop_feat. Not used as pseudogene evidence.
#   - premature_stop_codons column (from st:i / StopCodon=) = premature stops.
#   - missing_start_codon (da=-1) and clean_stop_codon (do=0) are new fields.
#   - mapq=255 means missing; mapq=0 is valid. Do not flag 0 as low quality.
#   - Donor/Acceptor on CDS = only written when NON-canonical.
#     has_noncanon_splice = TRUE means non-GT-AG detected.
#
# USAGE:
#   Rscript 03_report_pseudogenes.R \
#       --input  <outdir>/pseudogene_candidates.tsv \
#       --all    <outdir>/pseudogene_candidates_all_alignments.tsv \
#       --outdir <outdir>/report/
#
# DEPENDENCIES: tidyverse, optparse
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

# ==============================================================================
# 0. ARGUMENTS
# ==============================================================================

option_list <- list(
  make_option("--input",  type = "character"),
  make_option("--all",    type = "character", default = NULL),
  make_option("--outdir", type = "character", default = "report")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("--input is required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("=============================================================\n")
cat(" Pseudogene Candidate Report\n")
cat("=============================================================\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

pseudo <- read_tsv(opt$input, show_col_types = FALSE)
cat(sprintf("[03] Pseudogene candidates: %d\n\n", nrow(pseudo)))

all_aln <- NULL
if (!is.null(opt$all) && file.exists(opt$all)) {
  all_aln <- read_tsv(opt$all, show_col_types = FALSE)
  cat(sprintf("[03] Full alignment table:  %d rows\n\n", nrow(all_aln)))
}

# ==============================================================================
# 2. CLASSIFICATION SUMMARY
# ==============================================================================

label_map <- c(
  PROCESSED_PSEUDO    = "Processed pseudogene (retrotransposed)",
  NONPROCESSED_PSEUDO = "Non-processed pseudogene (gene duplicate)",
  UNITARY_PSEUDO      = "Unitary pseudogene (mutated in situ)",
  AMBIGUOUS           = "Ambiguous (insufficient evidence)"
)

cat("─────────────────────────────────────────────────────────────\n")
cat("SECTION 1 — Classification Counts\n")
cat("─────────────────────────────────────────────────────────────\n")

pseudo %>%
  count(classification) %>%
  arrange(desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         label = label_map[classification]) %>%
  { for (i in seq_len(nrow(.)))
      cat(sprintf("  %-52s  %5d  (%5.1f%%)\n",
                  .$label[i], .$n[i], .$pct[i])); . } %>%
  invisible()
cat("\n")

# ==============================================================================
# 3. PER-TYPE EVIDENCE SUMMARY
# ==============================================================================

cat("─────────────────────────────────────────────────────────────\n")
cat("SECTION 2 — Evidence Summary by Type\n")
cat("─────────────────────────────────────────────────────────────\n\n")

summarise_type <- function(df, label) {
  cat(sprintf("  ── %s (n = %d) ──\n", label, nrow(df)))
  if (nrow(df) == 0) { cat("     (none)\n\n"); return(invisible(NULL)) }

  cat(sprintf("    Identity:              mean=%.3f  median=%.3f  [%.3f–%.3f]\n",
    mean(df$identity, na.rm=TRUE), median(df$identity, na.rm=TRUE),
    min(df$identity,  na.rm=TRUE), max(df$identity, na.rm=TRUE)))

  cat(sprintf("    Query coverage:        mean=%.3f  median=%.3f\n",
    mean(df$q_coverage, na.rm=TRUE), median(df$q_coverage, na.rm=TRUE)))

  cat(sprintf("    Frameshifts:           total=%d  alignments with ≥1: %d (%.0f%%)\n",
    sum(df$frameshifts, na.rm=TRUE),
    sum(df$frameshifts > 0, na.rm=TRUE),
    100 * mean(df$frameshifts > 0, na.rm=TRUE)))

  # premature_stop_codons = st:i / StopCodon= (in-frame stops WITHIN alignment)
  cat(sprintf("    Premature stop codons: total=%d  alignments with ≥1: %d (%.0f%%)\n",
    sum(df$premature_stop_codons, na.rm=TRUE),
    sum(df$premature_stop_codons > 0, na.rm=TRUE),
    100 * mean(df$premature_stop_codons > 0, na.rm=TRUE)))

  # Missing start codon signal (da=-1)
  cat(sprintf("    Missing start codon:   %d (%.0f%%)  [da=-1]\n",
    sum(df$missing_start_codon, na.rm=TRUE),
    100 * mean(df$missing_start_codon, na.rm=TRUE)))

  # Clean stop codon = do=0 (proper termination)
  cat(sprintf("    Clean stop (do=0):     %d (%.0f%%)\n",
    sum(df$clean_stop_codon, na.rm=TRUE),
    100 * mean(df$clean_stop_codon, na.rm=TRUE)))

  # has_terminal_stop_feat = proper C-terminal termination (NOT pseudogene signal)
  cat(sprintf("    Terminal stop feature: %d (%.0f%%)  [C-terminal termination, good sign]\n",
    sum(df$has_terminal_stop_feat, na.rm=TRUE),
    100 * mean(df$has_terminal_stop_feat, na.rm=TRUE)))

  cat(sprintf("    Non-canonical splice:  %d (%.0f%%)\n",
    sum(df$has_noncanon_splice, na.rm=TRUE),
    100 * mean(df$has_noncanon_splice, na.rm=TRUE)))

  cat(sprintf("    CDS blocks:            mean=%.1f  single-exon (≤2): %d (%.0f%%)\n",
    mean(df$n_cds_blocks, na.rm=TRUE),
    sum(df$n_cds_blocks <= 2, na.rm=TRUE),
    100 * mean(df$n_cds_blocks <= 2, na.rm=TRUE)))

  cat(sprintf("    Span ratio:            mean=%.2f  compact (≤1.2): %d (%.0f%%)\n",
    mean(df$span_ratio, na.rm=TRUE),
    sum(df$span_ratio <= 1.2, na.rm=TRUE),
    100 * mean(df$span_ratio <= 1.2, na.rm=TRUE)))

  cat(sprintf("    Confidence:            mean=%.3f  high (≥0.7): %d (%.0f%%)\n",
    mean(df$confidence, na.rm=TRUE),
    sum(df$confidence >= 0.7, na.rm=TRUE),
    100 * mean(df$confidence >= 0.7, na.rm=TRUE)))
  cat("\n")
}

summarise_type(filter(pseudo, classification == "PROCESSED_PSEUDO"),    "PROCESSED PSEUDOGENES")
summarise_type(filter(pseudo, classification == "NONPROCESSED_PSEUDO"),  "NON-PROCESSED PSEUDOGENES")
summarise_type(filter(pseudo, classification == "UNITARY_PSEUDO"),       "UNITARY PSEUDOGENES")
summarise_type(filter(pseudo, classification == "AMBIGUOUS"),            "AMBIGUOUS")

# ==============================================================================
# 4. TOP CANDIDATES TABLE
# ==============================================================================

cat("─────────────────────────────────────────────────────────────\n")
cat("SECTION 3 — Top Candidates (confidence ≥ 0.7)\n")
cat("─────────────────────────────────────────────────────────────\n\n")

top <- pseudo %>%
  filter(classification != "AMBIGUOUS", confidence >= 0.7) %>%
  arrange(classification, desc(confidence)) %>%
  select(mrna_id, protein_id, contig, g_start, g_end, strand,
         classification, confidence, frameshifts, premature_stop_codons,
         missing_start_codon, clean_stop_codon,
         n_cds_blocks, has_noncanon_splice, span_ratio, identity, q_coverage)

if (nrow(top) > 0) {
  print(head(as.data.frame(top), 20), row.names = FALSE)
  cat(sprintf("\n  Showing %d of %d high-confidence candidates.\n\n",
              min(20, nrow(top)), nrow(top)))
} else {
  cat("  No candidates with confidence >= 0.7.\n\n")
}

# ==============================================================================
# 5. DIAGNOSTIC PLOTS
# ==============================================================================

cat("[03] Generating diagnostic plots...\n")

plot_file <- file.path(opt$outdir, "pseudogene_diagnostic_plots.pdf")
pdf(plot_file, width = 12, height = 8)

type_colours <- c(
  PROCESSED_PSEUDO    = "#E76F51",
  NONPROCESSED_PSEUDO = "#2A9D8F",
  UNITARY_PSEUDO      = "#E9C46A",
  AMBIGUOUS           = "#A8DADC",
  FUNCTIONAL          = "#457B9D",
  LOW_QUALITY_HIT     = "#BBBBBB"
)

# Helper: use all_aln if available, otherwise pseudo only
plot_df <- if (!is.null(all_aln)) all_aln else pseudo

# ── Plot 1: Classification count bar chart ────────────────────────────────────
p1 <- pseudo %>%
  count(classification) %>%
  mutate(classification = fct_reorder(classification, n)) %>%
  ggplot(aes(x = n, y = classification, fill = classification)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  scale_fill_manual(values = type_colours) +
  labs(title = "Pseudogene Classification Counts",
       subtitle = "Source: miniprot GFF alignment features",
       x = "Count", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank())
print(p1)

# ── Plot 2: Identity distribution by type ─────────────────────────────────────
p2 <- pseudo %>%
  filter(!is.na(identity)) %>%
  ggplot(aes(x = identity, fill = classification, colour = classification)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = type_colours) +
  scale_colour_manual(values = type_colours) +
  labs(title = "Alignment Identity Distribution by Pseudogene Type",
       subtitle = "Higher identity = more recently pseudogenised",
       x = "Identity (fraction of exact AA matches)", y = "Density",
       fill = "Type", colour = "Type") +
  theme_minimal(base_size = 13)
print(p2)

# ── Plot 3: Genomic span ratio vs identity ────────────────────────────────────
p3 <- pseudo %>%
  filter(!is.na(span_ratio), span_ratio < 20, !is.na(identity)) %>%
  ggplot(aes(x = span_ratio, y = identity,
             colour = classification, size = confidence)) +
  geom_point(alpha = 0.65) +
  geom_vline(xintercept = 1.2, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 1.25, y = 0.42,
           label = "No-intron threshold (1.2×)",
           hjust = 0, colour = "grey40", size = 3.2) +
  scale_colour_manual(values = type_colours) +
  scale_size_continuous(range = c(1.5, 5)) +
  labs(title = "Genomic Span Ratio vs. Alignment Identity",
       subtitle = "span_ratio = genomic_bp (incl. introns) / (protein_len × 3)\n~1.0 = no introns (processed pseudo); >>1.0 = introns retained",
       x = "Span ratio", y = "Identity",
       colour = "Type", size = "Confidence") +
  theme_minimal(base_size = 13)
print(p3)

# ── Plot 4: FRAMESHIFT count distribution ─────────────────────────────────────
p4 <- pseudo %>%
  mutate(fs_cap = factor(pmin(frameshifts, 5),
                         levels = 0:5, labels = c("0","1","2","3","4","5+"))) %>%
  count(classification, fs_cap) %>%
  ggplot(aes(x = fs_cap, y = n, fill = classification)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = type_colours) +
  labs(title = "Frameshift Count Distribution by Pseudogene Type",
       subtitle = "fs:i from PAF + Frameshift= from mRNA GFF attribute (max of two)",
       x = "Number of frameshifts", y = "Count", fill = "Type") +
  theme_minimal(base_size = 13)
print(p4)

# ── Plot 5: PREMATURE STOP CODON count distribution ──────────────────────────
# These are in-frame stop codons WITHIN the alignment (st:i / StopCodon=).
# Distinct from the stop_codon GFF feature which = proper C-terminal termination.
p5 <- pseudo %>%
  mutate(st_cap = factor(pmin(premature_stop_codons, 5),
                         levels = 0:5, labels = c("0","1","2","3","4","5+"))) %>%
  count(classification, st_cap) %>%
  ggplot(aes(x = st_cap, y = n, fill = classification)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = type_colours) +
  labs(title = "Premature In-Frame Stop Codon Count by Pseudogene Type",
       subtitle = "Source: st:i tag (PAF) and StopCodon= attribute (mRNA GFF)\nNOTE: this is NOT the stop_codon GFF feature, which signals proper C-terminal termination",
       x = "Number of premature in-frame stop codons", y = "Count", fill = "Type") +
  theme_minimal(base_size = 13)
print(p5)

# ── Plot 6: da/do ORF boundary signals ───────────────────────────────────────
# da = distance to nearest start codon (da=0 → ATG present; da=-1 → missing)
# do = distance to nearest stop codon  (do=0 → proper termination)
p6_data <- pseudo %>%
  mutate(
    start_status = case_when(
      missing_start_codon        ~ "Missing (da=-1)",
      as.numeric(da) == 0        ~ "At ATG (da=0)",
      TRUE                       ~ "Partial (da>0)"
    ),
    stop_status = case_when(
      clean_stop_codon           ~ "Clean (do=0)",
      TRUE                       ~ "Absent/far (do>0)"
    )
  )

p6a <- p6_data %>%
  count(classification, start_status) %>%
  ggplot(aes(x = classification, y = n, fill = start_status)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c(
    "Missing (da=-1)"  = "#E76F51",
    "At ATG (da=0)"    = "#2A9D8F",
    "Partial (da>0)"   = "#E9C46A"
  )) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Start Codon Status by Type (da tag)",
       subtitle = "da=-1: no upstream ATG found → disrupted start = pseudogene signal",
       x = NULL, y = "Fraction", fill = "Start status") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p6b <- p6_data %>%
  count(classification, stop_status) %>%
  ggplot(aes(x = classification, y = n, fill = stop_status)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c(
    "Clean (do=0)"        = "#2A9D8F",
    "Absent/far (do>0)"   = "#E76F51"
  )) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Stop Codon Status by Type (do tag)",
       subtitle = "do=0: alignment ends at stop codon → proper termination",
       x = NULL, y = "Fraction", fill = "Stop status") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Print as two side-by-side panels using cowplot if available, else separately
tryCatch({
  library(cowplot)
  print(plot_grid(p6a, p6b, ncol = 2,
                  labels = c("A: Start codon (da)", "B: Stop codon (do)")))
}, error = function(e) {
  print(p6a)
  print(p6b)
})

# ── Plot 7: Mutational profile heatmap ───────────────────────────────────────
# Shows the combination of disruption signals per type
p7 <- pseudo %>%
  mutate(
    has_frameshift   = frameshifts > 0,
    has_prem_stop    = premature_stop_codons > 0,
    has_miss_start   = missing_start_codon,
    has_no_clean_stop= !clean_stop_codon
  ) %>%
  pivot_longer(cols = c(has_frameshift, has_prem_stop,
                         has_miss_start, has_no_clean_stop),
               names_to = "signal", values_to = "present") %>%
  mutate(signal = recode(signal,
    has_frameshift    = "Frameshift (fs:i / Frameshift=)",
    has_prem_stop     = "Premature stop (st:i / StopCodon=)",
    has_miss_start    = "Missing start codon (da=-1)",
    has_no_clean_stop = "No clean stop (do>0)"
  )) %>%
  group_by(classification, signal) %>%
  summarise(pct_present = 100 * mean(present, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = signal, y = classification, fill = pct_present)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", pct_present)), size = 3.5) +
  scale_fill_gradient(low = "#f7f7f7", high = "#d62728",
                      name = "% with\nsignal") +
  labs(title = "Disruption Signal Prevalence by Pseudogene Type",
       subtitle = "Percentage of alignments showing each ORF disruption signal",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p7)

# ── Plot 8: CDS block count and splice site status ───────────────────────────
p8 <- pseudo %>%
  mutate(
    cds_cap = pmin(n_cds_blocks, 10),
    splice  = if_else(has_noncanon_splice, "Non-canonical splice", "Canonical/single-exon")
  ) %>%
  ggplot(aes(x = factor(cds_cap), fill = splice)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c(
    "Non-canonical splice"  = "#E9C46A",
    "Canonical/single-exon" = "#2A9D8F"
  )) +
  facet_wrap(~classification, scales = "free_y") +
  labs(title = "CDS Block Count and Splice Site Status",
       subtitle = "Non-canonical = Donor/Acceptor attr present on CDS (NN or non-GT-AG)\nCanonical GT-AG introns produce NO Donor/Acceptor attribute (per miniprot manual)",
       x = "CDS block count (capped at 10)", y = "Count", fill = "Splice sites") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(p8)

dev.off()
cat(sprintf("[03] Plots saved: %s  (8 panels)\n\n", plot_file))

# ==============================================================================
# 6. WRITE GFF3
# ==============================================================================

cat("[03] Writing pseudogene GFF3...\n")
gff3_file <- file.path(opt$outdir, "pseudogene_candidates.gff3")

so_type <- c(
  PROCESSED_PSEUDO    = "processed_pseudogene",
  NONPROCESSED_PSEUDO = "pseudogene",
  UNITARY_PSEUDO      = "pseudogene",
  AMBIGUOUS           = "pseudogene"
)

gff3_rows <- pseudo %>%
  filter(classification != "AMBIGUOUS") %>%
  mutate(
    gff_type     = so_type[classification],
    evidence_esc = str_replace_all(evidence, "[;\t]", " "),
    attributes   = paste0(
      "ID=",              mrna_id,                      ";",
      "Name=",            protein_id, "_pseudo;",
      "pseudo_type=",     classification,               ";",
      "confidence=",      confidence,                   ";",
      "identity=",        identity,                     ";",
      "q_coverage=",      q_coverage,                   ";",
      "frameshifts=",     frameshifts,                  ";",
      "premature_stops=", premature_stop_codons,         ";",
      "missing_start=",   missing_start_codon,          ";",
      "clean_stop=",      clean_stop_codon,             ";",
      "n_cds_blocks=",    n_cds_blocks,                 ";",
      "noncanon_splice=", has_noncanon_splice,          ";",
      "span_ratio=",      span_ratio,                   ";",
      "evidence=",        evidence_esc
    ),
    phase = "."
  ) %>%
  select(
    seqid      = contig,
    source     = classification,
    type       = gff_type,
    start      = g_start,
    end        = g_end,
    score      = confidence,
    strand     = strand,
    phase      = phase,
    attributes = attributes
  )

write_lines("##gff-version 3", gff3_file)
write_lines(paste0("# Pseudogene candidates — miniprot GFF — generated ",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
            gff3_file, append = TRUE)
write_tsv(gff3_rows, gff3_file, col_names = FALSE, append = TRUE)
cat(sprintf("[03] GFF3: %s  (%d features)\n\n", gff3_file, nrow(gff3_rows)))

# ==============================================================================
# 7. WRITE CLEAN TSV
# ==============================================================================

clean_file <- file.path(opt$outdir, "pseudogene_candidates_clean.tsv")
pseudo %>%
  filter(classification != "AMBIGUOUS", confidence >= 0.5) %>%
  arrange(classification, desc(confidence)) %>%
  write_tsv(clean_file)
cat(sprintf("[03] Clean TSV: %s\n", clean_file))

cat("\n=============================================================\n")
cat(" Report complete.\n")
cat("=============================================================\n")
