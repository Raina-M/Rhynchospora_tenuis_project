library(stats)

args <- commandArgs(trailingOnly = TRUE)
# mandatory input
AF_file <- args[1]
genomeFile <- args[2]
plotfile <- args[3]


########################### MAIN ###################
#
#
# read file
AF <- read.table(AF_file)
genome_sizes <- read.table(genomeFile)[,1:2]
chr_num <- dim(genome_sizes)[1]
max_chr <- max(genome_sizes[,2])

# prepare plots
pdf(file=plotfile, width=10, height=6)
par(mfrow = c(2, 1))

this_sample_co <- c()
for (chr in 1:chr_num){
  thisAF <- AF[AF$V1==chr,]
  chrsize <- genome_sizes[chr,2]
  
  # --------------- #
  # visualization   #
  # --------------- #
  
  # only for plotting
  ### plot 1 ###
  # AF with slight smoothing for a better visualization
  plot(thisAF$V2, thisAF$V4, type = "l", col = "blue", lwd = 1,
       main = paste("Chr", chr),
       xlab = "", ylab = "Normalized allele frequency",
       xlim=c(0, max_chr), ylim=c(-1,1),
       axes=F, cex.lab=1.5)
  axis(1, cex.axis=1.2,
       at=c(seq(0, chrsize, 2e+07), chrsize),
       labels = round(c(seq(0, chrsize, 2e+07), chrsize)/1e+06))
  axis(2, cex.axis=1.2, at = c(-1, 0, 1))
  #abline(h=0, lty=2, lwd=1.2)
  segments(x0=0,y0=0,x1=chrsize,y1=0, lty=2, lwd=1.2)
  mtext("Mb", at= chrsize+10e+06, padj = -1.35, side=1, line=2.0, cex=1.2)
}
dev.off()

