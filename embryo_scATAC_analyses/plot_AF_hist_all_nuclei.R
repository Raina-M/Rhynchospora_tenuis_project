WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/0_plots_for_MS_Fig_4g"

REC_em="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/ref_embryo_10XATAC/post_CO_analyses/allele_count_Cov_corrected_haps.txt"
REC_en="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/ref_endosperm_10XATAC/post_CO_analyses/allele_count_Cov_corrected_haps.txt"

PECP_em="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/PECP3_embryo_10XATAC/post_CO_analyses/allele_count_Cov.txt"
PECP_en="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/PECP3_endosperm_10XATAC/post_CO_analyses/allele_count_Cov.txt"

setwd(WD)

# read files - REC embryo and endosperm, PECP-48 embryo and endosperm
rec_em <- read.table(REC_em, header = F)
rec_en <- read.table(REC_en, header = F)
pecp_em <- read.table(PECP_em, header = F)
pecp_en <- read.table(PECP_en, header = F)

# merge the embryo and endosperm data together
# because em and en has contamination from each other
rec_chr1 <- c(rec_em$V2/(rec_em$V2 + rec_em$V3), rec_en$V2/(rec_en$V2+rec_en$V3))
rec_chr2 <- c(rec_em$V4/(rec_em$V4 + rec_em$V5), rec_en$V4/(rec_en$V4+rec_en$V5))
pecp_chr1 <- c(pecp_em$V2/(pecp_em$V2 + pecp_em$V3), pecp_en$V2/(pecp_en$V2+pecp_en$V3))
pecp_chr2 <- c(pecp_em$V4/(pecp_em$V4 + pecp_em$V5), pecp_en$V4/(pecp_en$V4+pecp_en$V5))

# plot the distribution of marker ratio: hap1/(hap1+hap2)
# fit to density
hist(rec_chr1, breaks = 100)
dens <- density(rec_chr1, from = 0, to = 1, n = 512)
lines(dens$x, dens$y*70, lwd = 1.2)
hist(rec_chr2, breaks = 100)

hist(pecp_chr1, breaks = 50)
hist(pecp_chr2, breaks = 50)

# A data array (one-dimension) contains the ratio data, so the value range is 0-1.
# Plot the distribution of this data and its density curve. Note that they are in the same plot, so scale the frequency of histogram to density curve. I expected two peaks/normal distributions. Decide the boundary of these two normal distribution using the local minimum between two peaks.

x <- pecp_chr2
# Compute kernel density
d <- density(x, from = 0, to = 1, n = 1000)

# Set pdf output
pdf("PECP_chr2_hist.pdf", width = 6, height = 4)

# Histogram scaled to density
hist(x, breaks = 50, freq = FALSE,
     xlim=c(0,1),
     col = "gray75",
     border = "white",
     xlab = "Ratio",
     main = "")

# Overlay density curve
lines(d, col = "#3c3c3c", lwd = 2)

# ---- Identify peaks and boundary ----
# Find local maxima (peaks)
is_peak <- function(y) diff(sign(diff(y))) == -2
peaks_idx <- which(is_peak(d$y)) + 1

# Only keep top 2 peaks (in case noise gives more)
top2 <- order(d$y[peaks_idx], decreasing = TRUE)[1:2]
peak_positions <- d$x[peaks_idx[top2]]
peak_positions <- sort(peak_positions)

# Find local minimum between the two peaks
range_idx <- which(d$x >= peak_positions[1] & d$x <= peak_positions[2])
min_idx <- range_idx[which.min(d$y[range_idx])]
boundary <- d$x[min_idx]

# Add peaks and boundary to plot
abline(v = peak_positions, col = c("#1b9e77","#1b9e77"), lwd = 2, lty = 2)
abline(v = boundary, col = "#d95f02", lwd = 2, lty = 3)

legend("topright",
       legend = c("Density", "Peaks", "Boundary"),
       col = c("#3c3c3c","#1b9e77","#d95f02"),
       lwd = c(2,2,2), lty = c(1,2,3))

# Add text beside lines
text(peak_positions[1], d$y[peaks_idx[top2]][1], 
     labels=sprintf("Peak1 = %.3f", peak_positions[1]),
     pos=4, col="#1b9e77", cex=0.9)

text(peak_positions[2], d$y[peaks_idx[top2]][2], 
     labels=sprintf("Peak2 = %.3f", peak_positions[2]),
     pos=4, col="#1b9e77", cex=0.9)

text(boundary, d$y[min_idx], 
     labels=sprintf("Boundary = %.3f", boundary),
     pos=4, col="#d95f02", cex=0.9)

cat("Peak positions:", peak_positions, "\n")
cat("Boundary (local minimum):", boundary, "\n")

dev.off()

