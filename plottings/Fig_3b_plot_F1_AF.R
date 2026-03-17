setwd("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/3_run_TIGER/")

# ------------------------------------------------------------
# Read marker allele count data
# Format: chrom pos hap1_allele hap1_count hap2_allele hap2_count
# ------------------------------------------------------------
df <- read.table("A03/input_corrected.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("chrom", "pos", "allele_h1", "count_h1", "allele_h2", "count_h2")

# Compute depth and allele frequency
df$depth <- df$count_h1 + df$count_h2
df <- df[df$depth>0,]
df$af_h1 <- df$count_h1/df$depth

# ------------------------------------------------------------
# Parameters for smoothing
# ------------------------------------------------------------
window_size <- 1000000    # 1 Mb window
step_size   <- 500000     # 500 kb step

# ------------------------------------------------------------
# Function to smooth AF by sliding window
# ------------------------------------------------------------
smooth_af <- function(chr_df, chrsize, window_size, step_size) {
  chr_df <- chr_df[order(chr_df$pos), ]
  starts <- seq(1, chrsize, by = step_size)
  
  smoothed <- data.frame(center = numeric(), mean_af = numeric())
  
  for (start in starts) {
    end <- start + window_size - 1
    idx <- which(chr_df$pos >= start & chr_df$pos <= end)
    if (length(idx) > 0) {
      smoothed <- rbind(smoothed,
                        data.frame(center = start + window_size / 2,
                                   mean_af = mean(chr_df$af_h1[idx], na.rm = TRUE)))
    }
  }
  return(smoothed)
}

# ------------------------------------------------------------
# Compute smoothed AF for each chromosome
# ------------------------------------------------------------
chrsizes <- read.table("/netscratch/dep_mercier/grp_marques/mzhang/HiC_maps/Rhynchospora_tenuis/references/Rtenuis.hap1.hifiasm.phased.hic.chrom.sizes")

chrom_list <- unique(df$chrom)
smoothed_list <- list()

for (ch in chrom_list) {
  dch <- df[df$chrom == ch, ]
  smoothed_list[[ch]] <- smooth_af(dch, chrsizes$V2[ch], window_size, step_size)
}

# ------------------------------------------------------------
# Set up plotting layout: main plot + right-side density
# ------------------------------------------------------------
pdf("A03_AF_test.pdf", width = 7.5, height=5)

layout(matrix(c(1, 2,
                3, 4), ncol = 2, byrow = TRUE),
       widths = c(4, 1.2), heights = c(1, 1))
par(oma = c(4, 4, 2, 1)) # outer margins

#cols <- c("#ff8010", "#1f77b4")
cols <- c("black", "black")

# ------------------------------------------------------------
# Loop over chromosomes for plotting
# ------------------------------------------------------------
for (i in seq_along(chrom_list)) {
  ch <- chrom_list[i]
  dch <- df[df$chrom == ch, ]
  sm <- smoothed_list[[ch]]
  chr_len <- chrsizes$V2[i]
  
  ## ---- Left panel: smoothed allele frequency ----
  par(mar = c(2, 4, 2, 0))
  plot(NA, NA, xlim = c(0, max(chrsizes$V2)), ylim = c(0, 1),
       xlab = "", ylab = "Allele frequency (hap1)",
       main = paste("Chr", ch), axes = FALSE)
  axis(1, at=c(seq(0, chr_len, 5e+07), chr_len),
       labels = round(c(seq(0, chr_len, 5e+07), chr_len)/1e+06), cex.axis=1.2)
  axis(2, at=seq(0,1,0.2), labels = c(0, "", "", "", "", 1))
  abline(h = c(0.25, 0.5, 0.75), lty = c(2, 2, 2), col = c("#3c3c3c", "#3c3c3c", "#3c3c3c"), lwd=1)
  lines(sm$center, sm$mean_af, col = cols[i], lwd = 1.5)
  
  ## ---- Right panel: rotated AF density ----
  par(mar = c(2, 0, 2, 4))
  dens <- density(dch$af_h1)
  x_density <- dens$y / max(dens$y) * chr_len * 0.15  # scale to same units as AF x-axis
  y_density <- dens$x
  
  plot(NA, NA, xlim = c(0, max(x_density)), ylim = c(0, 1),
       axes=F, xlab = "", ylab = "")
  lines(x_density, y_density, col = adjustcolor(cols[i], alpha.f = 0.75),
        border = cols[i], lwd = 1.5)
  axis(1, at=c(0,max(x_density)), labels = c(0, "Max"), las = 1)
  if (i == 2) mtext("AF density", side = 1, line = 2.5)
}

# ------------------------------------------------------------
# Add global x-axis label
# ------------------------------------------------------------
mtext("Position along chromosome (bp)", side = 1, outer = TRUE, line = 2.5)
dev.off()