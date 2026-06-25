
# Estimate divergence time based on mutation rate.

library(ggplot2)
library(scales)

setwd("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/compare_genomes/dN_dS_ratio/CoGe_Kn_Ks/estimate_divergence_time")
# ---------------------------
# User parameters
SNP_FILE <- "SNP_divergence_within_accession.txt"
KS_FILE <- "ALL_ACCESSION_Ks_purged.list"
LTR_FILE <- "ltr_insert_from_SubPhaser.txt"
BINS_Ks <- 50           # histogram bin count
BINS_SNP <- 50
BINS_LTR <- 50
MU <- 6.13e-9
# ---------------------------

######## Functions ###########

# Age conversion (Mya) helper
denom_note <- "t = d / (2·μ)"

div_to_mya <- function(d) (d / (2* MU)) / 1e6


# search local maximum based on density
local_maxima <- function(dens){
  # Input: output from density() 
  # Length encode the sign of difference
  rle <- rle(diff(as.vector(dens$y))>0)
  # starting point of searching local maximum
  starts <- cumsum(rle$lengths) - rle$lengths + 1
  # take the points where the rle is FALSE
  # i.e. difference goes from positive to negative
  maxima <- dens$x[starts[!rle$values]]
  return(maxima)  # x coordinates of the local maximum
}


########## End of Functions ###########


# Read Ks values
ks_raw <- read.table(KS_FILE, header = F)[,1]
ks <- log10(ks_raw)
ks <- ks[is.finite(ks)]

# read SNP data
snp <- read.table(SNP_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[,3]
snp <- snp[!is.na(snp)]/1e06

# read LTR data
ltr <- read.table(LTR_FILE, header = F)[,1]
ltr <- ltr[!is.na(ltr)]*1e06*2*1.3e-08  # LTR insertion time (million years) from SubPhaser uses mutation rate 1.3*1e08 

######################
# local maximal coordinates -> peaks of density plot -> potential important evolutionary events
ks_dens <- density(ks, kernel = "gaussian", adjust = 3)
ks_maxima <- local_maxima(ks_dens)

df_ks_dens <- data.frame(x=ks_dens$x, y=ks_dens$y)
ks_maxima_filtered <- df_ks_dens[df_ks_dens$x %in% ks_maxima & df_ks_dens$y>0.5, ]


################ Plot Ks ###############
# prepare dataframe for plot
df_ks <- data.frame(value=ks)

# plot
p_ks <- ggplot() +
  geom_histogram(aes(x = value, y = ..density..), data = df_ks,
                 bins=BINS_Ks, fill = alpha("steelblue", 0.25)) +
  geom_line(aes(x = x, y = y), data = df_ks_dens, color = "steelblue", linewidth=1) + 
  theme_bw() +
  xlab("log10(Ks)") +
  ylab("Density") +
  ggtitle("Ks distribution") +
  scale_x_continuous(limits=c(-4.3, 1.5), breaks=seq(-4,1,1)) +
  theme(plot.title = element_text(hjust = 0.5, size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        plot.margin = unit(c(1,1,1,1), "cm"))

# Annotate Ks component means and age ranges (convert means back to linear Ks)
for(i in seq_along(ks_maxima_filtered)) {
  xintercept <- ks_maxima_filtered$x[i]
  lab <- sprintf("~ %.1f Mya", div_to_mya(10^xintercept))
  p_ks <- p_ks +
    geom_vline(xintercept = xintercept, color = "#3E3E3E", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = xintercept, y = 0.85, label = lab, hjust = -0.2, angle = 0, size = 5, color = "#3E3E3E")
}

pdf("All_Acc_Ks_dens_plot.pdf",10,6)
print(p_ks)
dev.off()

################ Plot SNPs  ###############
# prepare dataframe for plot
df_snp <- data.frame(value=snp)

snp_dens <- density(snp, kernel = "gaussian", adjust = 3)
snp_maxima <- local_maxima(snp_dens)

df_snp_dens <- data.frame(x=snp_dens$x, y=snp_dens$y)
snp_maxima_filtered <- df_snp_dens[df_snp_dens$x %in% snp_maxima & df_snp_dens$y>0.5, ]

# plot
p_snp <- ggplot() +
  geom_histogram(aes(x = value, y = ..density..), data = df_snp,
                 bins=BINS_SNP, fill = alpha("firebrick", 0.25)) +
  geom_line(aes(x = x, y = y), data = df_snp_dens, color = "firebrick", linewidth=1) + 
  theme_bw() +
  xlab("# SNP / bp") +
  ylab("Density") +
  ggtitle("SNP density distribution") +
  #scale_x_continuous(limits=c(-4.3, 1.5), breaks=seq(-4,1,1)) +
  theme(plot.title = element_text(hjust = 0.5, size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        plot.margin = unit(c(1,1,1,1), "cm"))

# Annotate Ks component means and age ranges (convert means back to linear Ks)
for(i in seq_along(snp_maxima_filtered)) {
  xintercept <- snp_maxima_filtered$x[i]
  lab <- sprintf("~ %.1f Mya", div_to_mya(xintercept))
  p_snp <- p_snp +
    geom_vline(xintercept = xintercept, color = "#3E3E3E", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = xintercept, y = 500, label = lab, hjust = -0.2, angle = 0, size = 5, color = "#3E3E3E")
}

pdf("SNP_dens_plot.pdf",10,6)
print(p_snp)
dev.off()

################ Plot LTRs ###############
# prepare dataframe for plot
df_ltr <- data.frame(value=ltr)

ltr_dens <- density(ltr, kernel = "gaussian", adjust = 3)
ltr_maxima <- local_maxima(ltr_dens)

df_ltr_dens <- data.frame(x=ltr_dens$x, y=ltr_dens$y)
ltr_maxima_filtered <- df_ltr_dens[df_ltr_dens$x %in% ltr_maxima & df_ltr_dens$y>0.5, ]

# plot
p_ltr <- ggplot() +
  geom_histogram(aes(x = value, y = ..density..), data = df_ltr,
                 bins=BINS_LTR, fill = alpha("darkorange", 0.25)) +
  geom_line(aes(x = x, y = y), data = df_ltr_dens, color = "darkorange", linewidth=1) + 
  theme_bw() +
  xlab("LTR divergence") +
  ylab("Density") +
  ggtitle("LTR divergence distribution") +
  scale_x_continuous(limits=c(0, 0.2), breaks=seq(0,0.2,0.02)) +
  theme(plot.title = element_text(hjust = 0.5, size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        plot.margin = unit(c(1,1,1,1), "cm"))

# Annotate Ks component means and age ranges (convert means back to linear Ks)
for(i in seq_along(ltr_maxima_filtered)) {
  xintercept <- ltr_maxima_filtered$x[i]
  lab <- sprintf("~ %.1f Mya", div_to_mya(xintercept))
  p_ltr <- p_ltr +
    geom_vline(xintercept = xintercept, color = "#3E3E3E", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = xintercept, y = 50, label = lab, hjust = -0.2, angle = 0, size = 5, color = "#3E3E3E")
}

pdf("LTR_dens_plot.pdf",10,6)
print(p_ltr)
dev.off()




############## Plot: divergence time estimated by Ks vs. by SNP ns. by LTR ###############
# Calculate divergence time
h_ks <- hist(ks, breaks = BINS_Ks, plot = F)
df_ks_dens <- df_ks_dens[df_ks_dens$x>=min(h_ks$breaks) & df_ks_dens$x<=max(h_ks$breaks),]
df_ks_dens$age <- div_to_mya(10^(df_ks_dens$x))
df_ks_dens$scaled_density <- df_ks_dens$y/sum(df_ks_dens$y)
df_ks_dens <- df_ks_dens[df_ks_dens$age>0.05 & df_ks_dens$age<500,]

h_snp <- hist(snp, breaks = BINS_SNP, plot = F)
df_snp_dens <- df_snp_dens[df_snp_dens$x>=min(h_snp$breaks) & df_snp_dens$x<=max(h_snp$breaks),]
df_snp_dens$age <- div_to_mya(df_snp_dens$x)
df_snp_dens$scaled_density <- df_snp_dens$y/sum(df_snp_dens$y)
df_snp_dens <- df_snp_dens[df_snp_dens$age>0.05 & df_snp_dens$age<500,]

h_ltr <- hist(ltr, breaks = BINS_LTR, plot = F)
df_ltr_dens <- df_ltr_dens[df_ltr_dens$x>=min(h_ltr$breaks) & df_ltr_dens$x<=max(h_ltr$breaks),]
df_ltr_dens$age <- div_to_mya(df_ltr_dens$x)
df_ltr_dens$scaled_density <- df_ltr_dens$y/sum(df_ltr_dens$y)
df_ltr_dens <- df_ltr_dens[df_ltr_dens$age>0.05 & df_ltr_dens<500,]

# x interception
# Annotate age
mean_ks_age <- div_to_mya(10^(ks_maxima_filtered$x[1]))
mean_snp_age <- div_to_mya(snp_maxima_filtered$x[1])
mean_ltr_age <- div_to_mya(ltr_maxima_filtered$x[1])

lab_ks <- sprintf("~ %.1f Mya", mean_ks_age)
lab_snp <- sprintf("~ %.1f Mya", mean_snp_age)
lab_ltr <- sprintf("~ %.1f Mya", mean_ltr_age)

# plot
p <- ggplot() +
  geom_line(data = df_ks_dens, aes(x = age, y = scaled_density),
            color = "steelblue", linetype = "solid", linewidth = 1.1) +
  geom_line(data = df_snp_dens, aes(x = age, y = scaled_density),
            color = "firebrick", linetype = "solid", linewidth = 1.1) +
  geom_line(data = df_ltr_dens, aes(x = age, y = scaled_density),
            color = "darkorange", linetype = "solid", linewidth = 1.1) +
  geom_vline(xintercept = mean_ks_age, color = "steelblue", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = mean_snp_age, color = "firebrick", linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = mean_ltr_age, color = "darkorange", linetype = "dashed", size = 0.8) +
  annotate("text", x = mean_ks_age, y = 0.022, label = lab_ks, hjust = -0.1, angle = 0, size = 6, color = "steelblue") +
  annotate("text", x = mean_snp_age, y = 0.022, label = lab_snp, hjust = 1.1, angle = 0, size = 6, color = "firebrick") +
  annotate("text", x = mean_ltr_age, y = 0.022, label = lab_ltr, hjust = -0.1, angle = 0, size = 6, color = "darkorange") +
  theme_bw() +
  xlab("Time (Mya, log scale)") +
  ylab("Relative density") +
  ggtitle("Ks vs. LTR vs. SNP divergence") +
  scale_x_log10(breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(10^.x))) +
  annotation_logticks() +
  theme(plot.title = element_text(hjust = 0.5, size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        plot.margin = unit(c(1,1,1,1), "cm"))

pdf("All_acc_KS_SNP_LTR_dens_TIME.pdf",10,6)
print(p)
dev.off()


