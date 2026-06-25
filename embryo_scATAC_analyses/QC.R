setwd("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/ref_embryo_10XATAC/")

read_num <- read.table("read_num.stats", skip=1)[,1]
hist(read_num[read_num>=500 & read_num<3e+04], breaks = 50,
     main="",
     xlab = 'Number of reads in each cell',
     col = '#3C3C3C',
     border = 'black',
     axes = F,
     cex.lab=1.1)
axis(side = 1, col.axis = '#333333', col = NA, col.ticks = 1, cex.axis=1.1,
     at=seq(0, 50000, 5000))
axis(side = 2, col.axis = '#333333', col = NA, col.ticks = 1, cex.axis=1.1)

##### barcode rank plot #####
cutoff=10000
valid_reads <- sort(read_num[which(read_num>=cutoff)], decreasing = T)
bkgd_reads <- sort(read_num[which(read_num<cutoff)], decreasing = T)

# prepare output: update the path if necessary
pdf("barcode_rank_plot.pdf",
    family="Helvetica", height=6, width=6)
par(mai = c(1.2, 1, 1, 0.5)); # margin: bottom, left, top, right

# valid barcodes
plot(x=log10(seq(1, length(valid_reads), 1)),
     y=log10(sort(valid_reads, decreasing = T)),
     xlim = c(0,6), ylim=c(0,8),
     cex=0.5, col="blue",
     main = "Barcode Rank Plot",
     xlab = "Barcodes", ylab = 'Read counts',
     border = 'black', axes = F,
     cex.lab=1.5, cex.main=1.5)
# background barcodes
points(x=log10(seq(length(valid_reads)+1, length(read_num), 1)),
       y=log10(bkgd_reads),
       pch=1, cex=0.5, col="darkgrey")
box(lwd=2)
axis(side = 1, col.ticks = 1, cex.axis=1.25,
     at=seq(0,6,1), labels = c(1, 10, 100, 1000, "10k", "100k", "1M"))
axis(side = 2, col.ticks = 1, cex.axis=1.25,
     at=seq(0,8,2), labels = c(1, 100, "10k", "1M", "100M"))
legend(4, 7.5, legend=c("Cell", "Background"),
       col=c("blue", "darkgrey"), cex=1,
       box.lty=0, pch=1)
dev.off()



##### aligned read number and aln rates #####
df_aln <- read.table("aligned_num.stats")[,2:3]

# prepare output
pdf("alinged_reads_distr.pdf",
    family="Helvetica", height=11, width=6)
par(mai = c(0.8, 1, 0.5, 0.5)); # margin: bottom, left, top, right
m <- cbind(1:3)
layout(m, heights =rep(1,3))

# aln rates
hist(df_aln$V3/df_aln$V2, breaks = 50,
     main="Alignment rate distribution across cells",
     xlab="Alignment rate",
     ylab = 'Number of cells',
     col = '#3C3C3C',
     border = 'black',
     axes = F,
     cex.lab=1.1)
box(lwd=1.25)
axis(side = 1, col.ticks = 1, cex.axis=0.8,
     at=seq(0,1,0.1))
axis(side = 2, col.ticks = 1, cex.axis=0.8,
     at=seq(0,600,100))

# Number of uniquely aligned reads
hist(log10(df_aln[,1]), breaks = 50,
     main="Number of reads",
     xlab="# reads",
     ylab = 'Number of cells',
     col = '#3C3C3C',
     border = 'black',
     axes = F,
     cex.lab=1.1)
box(lwd=1.25)
axis(side = 1, col.ticks = 1, cex.axis=0.8,
     at=seq(4,8,1), labels=c('10k', '100k', '1M', '10M', '100M'))
axis(side = 2, col.ticks = 1, cex.axis=0.8,
     at=seq(0,25,5))
     
# Number of uniquely aligned reads
hist(log10(df_aln[,2]), breaks = 50,
     main="Number of aligned reads (unique mapping)",
     xlab="# aligned reads",
     ylab = 'Number of cells',
     col = '#3C3C3C',
     border = 'black',
     axes = F,
     cex.lab=1.1)
box(lwd=1.25)
axis(side = 1, col.ticks = 1, cex.axis=0.8,
     at=seq(0,6,1), labels=c(1, 10, 100, '1k', '10k', '100k', '1M'))
axis(side = 2, col.ticks = 1, cex.axis=0.8,
     at=seq(0,500,100))
#abline(v=log10(5000), col='blue', lty="dashed", lwd=1.25)

dev.off()


dim(df_aln[df_aln$V3/df_aln$V2 >= 0.65 & df_aln$V3>=10000,])


##### switch rates #####

# read file
sw <- read.table("switches.stats")[,2:3]

# prepare output
pdf("switch_rate.pdf",
    family="Helvetica", height=8, width=6)
par(mai = c(1, 1, 0.8, 0.5)); # margin: bottom, left, top, right
m <- cbind(1:2)
layout(m, heights =rep(1,2))

# marker number distribution
hist(sw[,1], breaks=50,
     xlab = "Marker number", ylab = "Count of cells",
     main="Marker number distribution",
     col = '#3C3C3C',
     border = 'black',
     cex.lab=1.1)
box(lwd=1.25)
# switch rate distribution
hist(sw$V3/sw$V2, breaks = 50,
     xlab = "Switch rate", ylab = "Count of cells",
     main="Switch rate distribution across cells",
     col = '#3C3C3C',
     border = 'black',
     axes = F,
     cex.lab=1.1)
box(lwd=1.25)
axis(side = 1, col.ticks = 1, cex.axis=0.8,
     at=seq(0,0.5,0.1))
axis(side = 2, col.ticks = 1, cex.axis=0.8,
     at=seq(0,200,50))

dev.off()     
     
# dot plot of switch rates
plot(x=sw$V2, y=sw$V3)
abline(a=0, b=0.07, col='red', lwd=1.2)
