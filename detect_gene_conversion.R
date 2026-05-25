################################################################################
# Gene Conversion Detection in Selfed F1 Plants 
# ------------------------------------------------------------------------------
# Key improvements over threshold-only approaches:
#   (A) Reference allele mapping bias correction
#         Illumina reads map preferentially to the reference allele, causing
#         the expected het frequency to be < 0.5.  We estimate p_het empirically
#         from putatively heterozygous markers and use it as the null for all
#         subsequent tests.
#   (B) Per-SNP two-sided binomial test
#         H0: hap2_count ~ Binomial(total_depth, p_het)
#         A significant result means this SNP's allele ratio is inconsistent
#         with heterozygosity at the bias-corrected expected fraction.
#   (C) Benjamini-Hochberg FDR correction
#         Applied globally across ALL testable SNPs.
#   (D) Merge adjacent significant SNPs into tracts (gap ≤ max_gap_bp)
#   (E) Depth-ratio + Wilcoxon test to distinguish gene conversion from
#         deletion / hemizygosity
#         Conversion preserves depth; deletion halves it (~0.5 × flanking).
#   (F) Local depth normalisation (GC-content / mappability bias)
#         A sliding-window median normalises depth before the depth-drop test,
#         removing the systematic low-depth artefact seen in high-GC windows.
#
# Other biases noted but requiring external data to correct fully:
#   * Strand bias  — needs per-allele F/R read counts (not in input).
#   * Paralogous sequence variants (PSVs) — reads from duplicated loci
#     mis-map here, mimicking het.  Filter by mapping quality upstream (MAPQ).
#   * Clustered sequencing errors in low-complexity / homopolymer regions —
#     pre-filter with a base-quality / RepeatMasker BED if available.
#   * Recombination / LOH tracts — long hom runs expected from selfing; they
#     are not gene conversions.  Enforcing max_gc_bp and max_gc_markers guards
#     against them being called as conversions.
#
# Input columns (whitespace-delimited, no header):
#   1 chrom  2 pos  3 hap1_allele  4 hap1_count  5 hap2_allele  6 hap2_count
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# ==============================================================================
# 1.  PARAMETERS
# ==============================================================================
PARAMS <- list(
  # Depth filter
  min_depth          = 5,       # Discard markers with total depth < this

  # Homo- and heterozygosity thresholds
  p_het_lo           = 0.25,
  p_het_hi           = 0.75,

  # Statistical testing
  fdr_alpha          = 0.05,    # Benjamini-Hochberg FDR threshold

  # Local-depth normalisation window (number of markers, must be odd)
  depth_window       = 51,

  # Tract merging
  max_gap_bp         = 5000,    # Merge sig. SNPs separated by ≤ this into one tract

  # Gene-conversion tract size limits
  max_gc_markers     = 20,      # Max markers within a converted tract
  max_gc_bp          = 5000,   # Max physical length (bp) of a converted tract
  min_gc_markers     = 2,       # Min markers within a converted tract
  min_gc_bp          = 5,       # Min physical length (bp) of a converted tract

  # Depth-drop test
  flank_n_markers    = 20,      # Flanking markers used for depth comparison
  depth_ratio_min    = 0.70,    # Depth ratio below this → suspect deletion
  depth_wilcox_alpha = 0.05     # Wilcoxon p threshold for depth-drop significance
)

# ==============================================================================
# 2.  DATA INPUT
# ==============================================================================
read_markers <- function(file_path) {
  read.table(file_path, header = FALSE, sep = "\t",
             col.names = c("chrom","pos","hap1_allele","hap1_count",
                           "hap2_allele","hap2_count"),
             colClasses = c("character","integer","character","integer",
                            "character","integer"))
}


# ==============================================================================
# 3.  REFERENCE ALLELE BIAS ESTIMATION
# ==============================================================================
# In Illumina alignments, hap2 (alt) allele reads map slightly less
# efficiently than hap1 (ref) reads, so the empirical het frequency is
# systematically below 0.5.  We estimate the true expected het fraction
# from markers whose raw hap2_freq falls within a broad het window.
# Returns: scalar p_het (to replace 0.5 in all downstream tests).

estimate_p_het <- function(df,
                           min_depth = PARAMS$min_depth,
                           lo = PARAMS$p_het_lo,
                           hi = PARAMS$p_het_hi) {
  cands <- df %>%
    filter(total_depth >= min_depth,
           hap2_freq  >= lo,
           hap2_freq  <= hi)

  if (nrow(cands) < 30) {
    warning("< 30 candidate het markers; using p_het = 0.5.")
    return(0.5)
  }

  p_hat <- median(cands$hap2_freq, na.rm = TRUE)

  cat(sprintf(paste0(
    "--- Reference allele bias ---\n",
    "  Candidate het markers : %d\n",
    "  Estimated p_het       : %.4f  (bias from 0.5: %+.4f)\n\n"),
    nrow(cands), p_hat, p_hat - 0.5))

  p_hat
}

# ==============================================================================
# 4.  PER-SNP BINOMIAL TEST  +  BH FDR CORRECTION
# ==============================================================================
# Fast vectorised two-sided binomial p-value (doubles the smaller tail).
binom2_pval <- function(x, n, p) {
  pmin(2 * pmin(pbinom(x, n, p), pbinom(x - 1L, n, p, lower.tail = FALSE)), 1)
}

# Annotates df with:
#   p_raw         raw binomial p-value (NA if depth < min_depth)
#   q_bh          BH-corrected q-value
#   sig_deviation logical: q_bh < fdr_alpha
#   deviation_dir "hap1_excess" | "hap2_excess" | "none"

add_binom_test <- function(df, p_het,
                            min_depth = PARAMS$min_depth,
                            fdr_alpha = PARAMS$fdr_alpha) {
  df <- df %>%
    mutate(
      p_raw = ifelse(
        total_depth >= min_depth,
        binom2_pval(hap2_count, total_depth, p_het),
        NA_real_)
    )

  # BH correction is applied globally across ALL chromosomes at once.
  testable <- which(!is.na(df$p_raw))
  q_vals   <- rep(NA_real_, nrow(df))
  q_vals[testable] <- p.adjust(df$p_raw[testable], method = "BH")
  df$q_bh  <- q_vals

  df <- df %>%
    mutate(
      sig_deviation = !is.na(q_bh) & q_bh < fdr_alpha,
      deviation_dir = case_when(
        !sig_deviation         ~ "none",
        hap2_freq <  p_het     ~ "hap1_excess",
        hap2_freq >  p_het     ~ "hap2_excess",
        TRUE                   ~ "none"
      )
    )
  df
}

# ==============================================================================
# 5.  LOCAL DEPTH NORMALISATION  (GC-content / mappability bias correction)
# ==============================================================================
# A sliding-window median of total depth captures local sequencing
# non-uniformity driven by GC content, repetitive elements, and mappability.
# Dividing each marker's depth by this local baseline gives depth_norm ≈ 1.0
# everywhere in diploid regions and ≈ 0.5 in hemizygous regions, regardless
# of raw depth variation across the genome.

local_depth_normalize <- function(df, window = PARAMS$depth_window) {
  half <- (window - 1L) %/% 2L

  df %>%
    group_by(chrom) %>%
    arrange(pos, .by_group = TRUE) %>%
    mutate(
      local_med_depth = {
        n   <- n()
        dep <- total_depth
        vapply(seq_len(n), function(i) {
          idx <- max(1L, i - half):min(n, i + half)
          median(dep[idx], na.rm = TRUE)
        }, numeric(1))
      },
      depth_norm = ifelse(local_med_depth > 0,
                          total_depth / local_med_depth, NA_real_)
    ) %>%
    ungroup()
}

# ==============================================================================
# 6.  MERGE SIGNIFICANT SNPs INTO CANDIDATE TRACTS
# ==============================================================================
# Adjacent significant SNPs separated by ≤ max_gap_bp are merged into one
# tract.  The resulting tracts are then filtered for size (gene-conversion
# tracts are short).

merge_sig_snps <- function(df, max_gap_bp = PARAMS$max_gap_bp) {
  df %>%
    filter(sig_deviation) %>%
    group_by(chrom) %>%
    arrange(pos, .by_group = TRUE) %>%
    mutate(
      gap      = pos - lag(pos, default = first(pos)),
      new_grp  = row_number() == 1L | gap > max_gap_bp,
      tract_id = paste0(chrom, "_T", formatC(cumsum(new_grp), width = 4,
                                              flag = "0"))
    ) %>%
    ungroup()
}

# ==============================================================================
# 7.  BUILD TRACT TABLE + FLANKING CONTEXT CHECK
# ==============================================================================
# Condenses each merged tract to one summary row; keeps only short tracts
# whose internal deviation direction differs from at least one flanking block
# (the hallmark of a localised conversion event vs a long LOH block).

get_flanking_state <- function(chr, s_pos, e_pos, df_full, flank_n) {
  chr_df   <- df_full %>% filter(chrom == chr) %>% arrange(pos)
  left_df  <- chr_df  %>% filter(pos < s_pos)  %>% tail(flank_n)
  right_df <- chr_df  %>% filter(pos > e_pos)  %>% head(flank_n)

  dom_dir <- function(x) {
    x <- x %>% filter(!is.na(deviation_dir))
    if (nrow(x) == 0L) return("none")
    names(sort(table(x$deviation_dir), decreasing = TRUE))[1]
  }

  list(
    left_dom_dir     = dom_dir(left_df),
    right_dom_dir    = dom_dir(right_df),
    left_mean_depth  = mean(left_df$total_depth,  na.rm = TRUE),
    right_mean_depth = mean(right_df$total_depth, na.rm = TRUE),
    left_mean_freq   = mean(left_df$hap2_freq,    na.rm = TRUE),
    right_mean_freq  = mean(right_df$hap2_freq,   na.rm = TRUE)
  )
}

build_tract_table <- function(sig_snps, df_full,
                               max_gc_markers = PARAMS$max_gc_markers,
                               max_gc_bp      = PARAMS$max_gc_bp,
                               min_gc_markers = PARAMS$min_gc_markers,
                               min_gc_bp      = PARAMS$min_gc_bp,
                               flank_n        = PARAMS$flank_n_markers) {
  if (nrow(sig_snps) == 0L) return(tibble())

  tracts <- sig_snps %>%
    group_by(chrom, tract_id) %>%
    summarise(
      start_pos       = min(pos),
      end_pos         = max(pos),
      n_markers       = n(),
      tract_bp        = end_pos - start_pos,
      mean_hap2_freq  = mean(hap2_freq,    na.rm = TRUE),
      mean_depth      = mean(total_depth,  na.rm = TRUE),
      mean_depth_norm = mean(depth_norm,   na.rm = TRUE),
      mean_q_bh       = mean(q_bh,         na.rm = TRUE),
      dominant_dir    = {
        tb <- sort(table(deviation_dir), decreasing = TRUE)
        names(tb)[1]
      },
      .groups = "drop"
    ) %>%
    filter(n_markers <= max_gc_markers, tract_bp <= max_gc_bp,
           n_markers >= min_gc_markers, tract_bp >= min_gc_bp)

  if (nrow(tracts) == 0L) return(tibble())

  # Flanking context (rowwise — one API call per tract)
  tracts <- tracts %>%
    rowwise() %>%
    mutate(
      .fl = list(get_flanking_state(chrom, start_pos, end_pos, df_full, flank_n))
    ) %>%
    unnest_wider(.fl) %>%
    ungroup()

  # Keep only tracts where internal direction ≠ at least one flanking direction
  tracts %>%
    filter(left_dom_dir != dominant_dir | right_dom_dir != dominant_dir) %>%
    mutate(
      # Confidence: BOTH flanks same and different from tract = high confidence
      both_flanks_agree = (left_dom_dir == right_dom_dir) &
                          (left_dom_dir != dominant_dir),
      gc_confidence = ifelse(both_flanks_agree, "HIGH", "LOW")
    )
}

# ==============================================================================
# 8.  DEPTH TEST: GENE CONVERSION  vs.  DELETION / HEMIZYGOSITY
# ==============================================================================
# For each candidate tract we compare its total read depth to the flanking
# regions with two complementary tests:
#   (i)  depth_ratio = mean(tract_depth) / mean(flanking_depth)
#           Ratio ≈ 1.0 → diploid (conversion);  ≈ 0.5 → hemizygous deletion
#   (ii) Wilcoxon rank-sum test (one-sided: H1: tract < flanking)
#           Significant → depth is genuinely lower in the tract
# A tract is flagged as deletion / hemizygosity if EITHER test fires.

run_depth_test <- function(chr, s_pos, e_pos, df_full, flank_n) {
  chr_df <- df_full %>% filter(chrom == chr, total_depth > 0) %>% arrange(pos)

  tract_d   <- chr_df %>% filter(pos >= s_pos, pos <= e_pos) %>%
    pull(total_depth)
  left_d    <- chr_df %>% filter(pos < s_pos)  %>% tail(flank_n) %>%
    pull(total_depth)
  right_d   <- chr_df %>% filter(pos > e_pos)  %>% head(flank_n) %>%
    pull(total_depth)
  flank_d   <- c(left_d, right_d)

  if (length(tract_d) == 0L || length(flank_d) < 2L)
    return(list(depth_ratio = NA_real_, wilcox_p = NA_real_))

  list(
    depth_ratio = mean(tract_d) / mean(flank_d),
    wilcox_p    = tryCatch(
      wilcox.test(tract_d, flank_d, alternative = "less")$p.value,
      error = function(e) NA_real_)
  )
}

depth_test_tracts <- function(tracts, df_full,
                               depth_ratio_min    = PARAMS$depth_ratio_min,
                               depth_wilcox_alpha = PARAMS$depth_wilcox_alpha,
                               flank_n            = PARAMS$flank_n_markers) {
  if (nrow(tracts) == 0L) return(tracts)

  tracts %>%
    rowwise() %>%
    mutate(.dt = list(run_depth_test(chrom, start_pos, end_pos, df_full, flank_n))) %>%
    unnest_wider(.dt) %>%
    ungroup() %>%
    mutate(
      depth_suspect = (!is.na(depth_ratio) & depth_ratio < depth_ratio_min) |
                      (!is.na(wilcox_p)    & wilcox_p    < depth_wilcox_alpha),
      call = ifelse(depth_suspect,
                    "DELETION/HEMIZYGOSITY",
                    "GENE CONVERSION")
    )
}

# ==============================================================================
# 9.  REPORTING
# ==============================================================================
report_results <- function(tracts, df_full, p_het) {
  n_tested  <- sum(!is.na(df_full$p_raw))
  n_sig     <- sum(df_full$sig_deviation, na.rm = TRUE)
  n_tracts  <- nrow(tracts)
  n_gc      <- if (n_tracts > 0) sum(tracts$call == "GENE CONVERSION",      na.rm=TRUE) else 0L
  n_del     <- if (n_tracts > 0) sum(tracts$call == "DELETION/HEMIZYGOSITY", na.rm=TRUE) else 0L

  cat(paste0(
    "====================================================\n",
    "   Gene Conversion Detection — Statistical Report\n",
    "====================================================\n",
    sprintf("  Total markers            : %d\n",      nrow(df_full)),
    sprintf("  Markers tested           : %d\n",      n_tested),
    sprintf("  FDR-significant markers  : %d\n",      n_sig),
    sprintf("  Estimated p_het          : %.4f\n",    p_het),
    sprintf("  Candidate tracts (short) : %d\n",      n_tracts),
    sprintf("    → Gene conversions     : %d\n",      n_gc),
    sprintf("    → Deletions/hemizygous : %d\n\n",    n_del)
  ))

  if (n_tracts == 0L) { cat("No candidate tracts detected.\n"); return(invisible(NULL)) }

  fmt <- function(t, label) {
    cat(sprintf("--- %s ---\n", label))
    print(t %>% select(chrom, start_pos, end_pos, tract_bp, n_markers,
                        dominant_dir, mean_hap2_freq,
                        depth_ratio, wilcox_p,
                        gc_confidence, call) %>%
            mutate(across(c(start_pos, end_pos, tract_bp),
                           ~ format(.x, big.mark = ",")),
                   across(c(mean_hap2_freq, depth_ratio, wilcox_p),
                           ~ round(.x, 4))),
          n = Inf, width = 120)
    cat("\n")
  }

  gc_df  <- tracts %>% filter(call == "GENE CONVERSION")
  del_df <- tracts %>% filter(call == "DELETION/HEMIZYGOSITY")
  if (nrow(gc_df)  > 0L) fmt(gc_df,  "Gene Conversion Events")
  if (nrow(del_df) > 0L) fmt(del_df, "Deletion / Hemizygosity Events")

  invisible(tracts)
}

# ==============================================================================
# 10.  VISUALISATION
# ==============================================================================

# --- Plot A: allele frequency coloured by statistical significance --------
plot_allele_freq <- function(df, tracts = NULL, p_het = 0.5, chroms = NULL) {
  pd <- df %>%
    filter(!is.na(hap2_freq), total_depth >= PARAMS$min_depth) %>%
    { if (!is.null(chroms)) filter(., chrom %in% chroms) else . } %>%
    mutate(pt_class = case_when(
      sig_deviation & deviation_dir == "hap1_excess" ~ "Sig: hap1 excess",
      sig_deviation & deviation_dir == "hap2_excess" ~ "Sig: hap2 excess",
      TRUE ~ "Non-significant"
    )) %>%
    # CRITICAL FIX: Sorts FALSE (0) to the top of the dataframe, TRUE (1) to the bottom.
    # This guarantees significant points are plotted LAST, keeping them on top visually.
    arrange(sig_deviation)
  
  cols <- c("Non-significant"  = "grey72",
            "Sig: hap1 excess" = "#ff8010",
            "Sig: hap2 excess" = "#1f77b4")
  
  sizes <- c("Non-significant"  = 0.5,
             "Sig: hap1 excess" = 0.8,
             "Sig: hap2 excess" = 0.8)
  
  p <- ggplot(pd, aes(pos / 1e6, hap2_freq, colour = pt_class)) +
    geom_hline(yintercept = p_het, linetype = "dashed",
               colour = "black", linewidth = 0.55) +
    geom_point(alpha = 0.75) +
    scale_colour_manual(values = cols, name = NULL,
                        breaks = c("Non-significant", "Sig: hap1 excess", "Sig: hap2 excess")) +
    scale_size_manual(values = sizes, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    # MODIFIED: ncol = 1 (one per row), default fixed scale maps physical lengths accurately
    facet_wrap(~ chrom, ncol = 1,
               labeller = labeller(chrom = ~ paste0("Chr ", .x))) +
    labs(
      title    = "Hap2 allele frequency — statistical view",
      subtitle = sprintf(
        "Black dashed = p_het (%.3f) | Orange = hap1-excess | Blue = hap2-excess (FDR < %.2f)",
        p_het, PARAMS$fdr_alpha),
      x = "Position (Mb)", y = "Hap2 allele frequency"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"))
  
  p <- .add_tract_rects(p, tracts, chroms)
  p
}

# --- Plot B: locally-normalised depth ----------------------------------------
plot_depth_norm <- function(df, tracts = NULL, chroms = NULL) {
  pd <- df %>%
    filter(!is.na(depth_norm)) %>%
    { if (!is.null(chroms)) filter(., chrom %in% chroms) else . }
  
  p <- ggplot(pd, aes(pos / 1e6, depth_norm)) +
    geom_point(alpha = 0.40, size = 0.9, colour = "steelblue") +
    geom_hline(yintercept = 1.0, colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, colour = "darkorange",
               linetype = "dashed", linewidth = 0.5) +
    # MODIFIED: ncol = 1 (one per row)
    facet_wrap(~ chrom, ncol = 1,
               labeller = labeller(chrom = ~ paste0("Chr ", .x))) +
    labs(
      title    = "Locally normalised read depth",
      subtitle = "Orange dashed = 0.5× (hemizygous). Depth drop alongside allele shift = deletion.",
      x = "Position (Mb)", y = "Normalised depth"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"))
  
  p <- .add_tract_rects(p, tracts, chroms)
  p
}

# Internal helper: overlay tract rectangles + labels
.add_tract_rects <- function(p, tracts, chroms) {
  if (is.null(tracts) || nrow(tracts) == 0L) return(p)
  td <- tracts
  if (!is.null(chroms)) td <- td %>% filter(chrom %in% chroms)
  if (nrow(td) == 0L) return(p)
  
  # Calculate midpoint position in Mb for precise placement
  td <- td %>% mutate(mid_pos = (start_pos + end_pos) / 2 / 1e6)
  
  # Isolate events to apply distinct, hardcoded styling safely
  gc_df  <- td %>% filter(call == "GENE CONVERSION")
  del_df <- td %>% filter(call == "DELETION/HEMIZYGOSITY")
  
  # --- Overlay Gene Conversion Events (Bright Bold Red) ---
  if (nrow(gc_df) > 0L) {
    p <- p +
      # Bright red vertical line intersecting the entire Y-axis
      geom_vline(data = gc_df, aes(xintercept = mid_pos),
                 color = "red", linewidth = 0.7, alpha = 0.5) #+
      # High-visibility solid red badge at the top of the plot
      #geom_label(data = gc_df, aes(x = mid_pos, y = Inf, 
      #                             label = paste0("GC: ", gc_confidence)),
      #           vjust = 1.2, size = 3.8, fontface = "bold", 
      #           color = "white", fill = "red", inherit.aes = FALSE)
  }
  
  # --- Overlay Deletion Events (Orange-Red) ---
  if (nrow(del_df) > 0L) {
    p <- p +
      # Orange-red vertical line intersecting the entire Y-axis
      geom_vline(data = del_df, aes(xintercept = mid_pos),
                 color = "purple", linewidth = 0.8, alpha = 0.75) +
      # High-visibility solid orange-red badge at the top of the plot
      geom_label(data = del_df, aes(x = mid_pos, y = Inf, label = "DEL"),
                 vjust = 1.2, size = 3.8, fontface = "bold", 
                 color = "white", fill = "purple", inherit.aes = FALSE)
  }
  
  p
}
# ==============================================================================
# 11.  MAIN PIPELINE
# ==============================================================================

## ---- Load data ---------------------------------------------------------------
# read input arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
outdir <- args[2]

markers <- read_markers(paste0("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/3_run_TIGER/", sample, "/input_corrected.txt")) %>%
  mutate(total_depth = hap1_count + hap2_count,
         hap2_freq   = ifelse(total_depth > 0,
                              hap2_count / total_depth, NA_real_))

## ---- Step 1: Reference allele bias estimate ----------------------------------
p_het <- estimate_p_het(markers)

## ---- Step 2: Per-SNP binomial test + BH FDR correction ----------------------
markers <- add_binom_test(markers, p_het)

## ---- Step 3: Local depth normalisation (GC / mappability bias) ---------------
markers <- local_depth_normalize(markers)

## ---- Step 4: Merge significant SNPs → candidate tracts ----------------------
sig_snps   <- merge_sig_snps(markers)
tract_raw  <- build_tract_table(sig_snps, markers)

## ---- Step 5: Depth test — conversion vs deletion / hemizygosity -------------
tracts <- depth_test_tracts(tract_raw, markers)

## ---- Report ------------------------------------------------------------------
report_results(tracts, markers, p_het)

## ---- Plots -------------------------------------------------------------------
n_chr <- n_distinct(markers$chrom)
h     <- max(4, 4 * ceiling(n_chr / 2))

p_af <- plot_allele_freq(markers, tracts, p_het)
p_dp <- plot_depth_norm(markers, tracts)
#print(p_af)
#print(p_dp)

ggsave(paste0(outdir, "/", sample, "_gc_allele_freq.png"),  p_af, width = 14, height = h, dpi=300)
ggsave(paste0(outdir, "/", sample, "_gc_depth_norm.png"),   p_dp, width = 14, height = h, dpi=300)

## ---- Export ------------------------------------------------------------------
out_markers <- markers %>%
  select(chrom, pos, hap1_allele, hap1_count, hap2_allele, hap2_count,
         total_depth, hap2_freq, local_med_depth, depth_norm,
         p_raw, q_bh, sig_deviation, deviation_dir)

write.table(out_markers, paste0(outdir, "/", sample, "_markers_annotated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(tracts) && nrow(tracts) > 0L) {
  write.table(tracts, paste0(outdir, "/", sample, "_gc_tracts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Tracts written to ", paste0(outdir, "/", sample, "_gc_tracts.tsv"), "\n")
}

cat("Annotated markers written to ", paste0(outdir, "/", sample, "_markers_annotated.tsv"), "\n")
cat("Done.\n")
