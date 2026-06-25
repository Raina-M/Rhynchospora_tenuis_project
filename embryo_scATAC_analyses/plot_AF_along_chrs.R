setwd("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/0_plots_for_MS_Fig_4g/")

source="PECP_EN"  # REC embryo
# REC chr sizes
#chrsizes <- read.table("/netscratch/dep_mercier/grp_marques/mzhang/HiC_maps/Rhynchospora_tenuis/references/Rtenuis.hap1.hifiasm.phased.hic.chrom.sizes")
# PECP-48 chr sizes
chrsizes <- read.table("/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344C.PECP3.hap1.chr.fasta.fai")

# ------------------------------------------------------------
# Set up plotting layout: main plot + right-side density
# ------------------------------------------------------------
pdf(paste("AF_along_chrs_", source, ".pdf", sep=""), width = 6, height=5)

layout(matrix(c(1, 2), nrow = 2, byrow = TRUE),
       widths = c(6, 6), heights = c(1.2, 1.2))
par(oma = c(1, 1, 1, 1)) # outer margins

#cols <- c("#ff8010", "#1f77b4")
cols <- c("black", "black")

# ------------------------------------------------------------
# Loop over chromosomes for plotting
# ------------------------------------------------------------
for (chr in 1:2) {
  # read AF of this chromosome
  df <- read.table(paste(source, "_chr", chr, "_h1_AF.txt", sep = ""),
                   header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("pos", "af_h1")
  
  chr_len <- chrsizes$V2[chr]
  
  par(mar = c(2, 4, 2, 0))
  smoothingSpline = smooth.spline(df$pos, df$af_h1, spar=0.5)
  plot(NA, NA, xlim = c(0, max(chrsizes$V2)), ylim = c(0, 1),
       xlab = "", ylab = "Allele frequency (hap1)",
       main = paste("Chr", chr), axes = FALSE)
  axis(1, at=c(seq(0, chr_len, 5e+07), chr_len),
       labels = round(c(seq(0, chr_len, 5e+07), chr_len)/1e+06))
  axis(2, at=seq(0,1,0.5))
  abline(h = c(1/3, 0.5, 2/3), lty = c(2, 2, 2), col = c("#3c3c3c", "#3c3c3c", "#3c3c3c"), lwd=1)
  lines(smoothingSpline, lwd=2)
}

# ------------------------------------------------------------
# Add global x-axis label
# ------------------------------------------------------------
mtext("Position along chromosome (bp)", side = 1, outer = TRUE, line = 2.5)
dev.off()