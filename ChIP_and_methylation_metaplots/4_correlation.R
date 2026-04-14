
# Load libraries
library(data.table)
library(ggplot2)

setwd("/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/DNA_methylation/compare_with_pubera_breviuscula")
#  Load column mean of matrix
pubera_cenh3 <- t(read.table("pubera_ChIP/CENH3_Tyba_colMeans.txt"))
pubera_h3k4me3 <- t(read.table("pubera_ChIP/H3K4me3_Tyba_colMeans.txt"))
pubera_h3k9me2 <- t(read.table("pubera_ChIP/H3K9me2_Tyba_colMeans.txt"))
pubera_h3k27me3 <- t(read.table("pubera_ChIP/H3K27me3_Tyba_colMeans.txt"))
pubera_cpg <- t(read.table("pubera_DNA_methylation/CpG_Tyba_colMeans.txt"))
pubera_chg <- t(read.table("pubera_DNA_methylation/CHG_Tyba_colMeans.txt"))
pubera_chh <- t(read.table("pubera_DNA_methylation/CHH_Tyba_colMeans.txt"))

brevi_cenh3 <- t(read.table("breviuscula_ChIP/CENH3_Tyba_colMeans.txt"))
brevi_h3k4me3 <- t(read.table("breviuscula_ChIP/H3K4me3_Tyba_colMeans.txt"))
brevi_h3k9me2 <- t(read.table("breviuscula_ChIP/H3K9me2_Tyba_colMeans.txt"))
brevi_h3k27me3 <- t(read.table("breviuscula_ChIP/H3K27me3_Tyba_colMeans.txt"))
brevi_cpg <- t(read.table("breviuscula_DNA_methylation/CpG_Tyba_colMeans.txt"))
brevi_chg <- t(read.table("breviuscula_DNA_methylation/CHG_Tyba_colMeans.txt"))
brevi_chh <- t(read.table("breviuscula_DNA_methylation/CHH_Tyba_colMeans.txt"))

tenuis_cenh3 <- t(read.table("tenuis_ChIP/CENH3_Tyba_colMeans.txt"))
tenuis_h3k4me3 <- t(read.table("tenuis_ChIP/H3K4me3_Tyba_colMeans.txt"))
tenuis_h3k9me2 <- t(read.table("tenuis_ChIP/H3K9me2_Tyba_colMeans.txt"))
tenuis_h3k27me3 <- t(read.table("tenuis_ChIP/H3K27me3_Tyba_colMeans.txt"))
tenuis_cpg <- t(read.table("tenuis_DNA_methylation/CpG_Tyba_colMeans.txt"))
tenuis_chg <- t(read.table("tenuis_DNA_methylation/CHG_Tyba_colMeans.txt"))
tenuis_chh <- t(read.table("tenuis_DNA_methylation/CHH_Tyba_colMeans.txt"))


# --- TEST 1: Correlation (Pattern/Shape Similarity) ---
# CENH3
cor_pb_cenh3 <- cor(pubera_cenh3, brevi_cenh3, method = "spearman")
cor_pt_cenh3 <- cor(pubera_cenh3, tenuis_cenh3, method = "spearman")
cor_bt_cenh3 <- cor(brevi_cenh3, tenuis_cenh3, method = "spearman")

# H3K4me3
cor_pb_h3k4me3 <- cor(pubera_h3k4me3, brevi_h3k4me3, method = "spearman")
cor_pt_h3k4me3 <- cor(pubera_h3k4me3, tenuis_h3k4me3, method = "spearman")
cor_bt_h3k4me3 <- cor(brevi_h3k4me3, tenuis_h3k4me3, method = "spearman")

# H3K9me2
cor_pb_h3k9me2 <- cor(pubera_h3k9me2, brevi_h3k9me2, method = "spearman")
cor_pt_h3k9me2 <- cor(pubera_h3k9me2, tenuis_h3k9me2, method = "spearman")
cor_bt_h3k9me2 <- cor(brevi_h3k9me2, tenuis_h3k9me2, method = "spearman")

# H3K27me3
cor_pb_h3k27me3 <- cor(pubera_h3k27me3, brevi_h3k27me3, method = "spearman")
cor_pt_h3k27me3 <- cor(pubera_h3k27me3, tenuis_h3k27me3, method = "spearman")
cor_bt_h3k27me3 <- cor(brevi_h3k27me3, tenuis_h3k27me3, method = "spearman")

# CpG
cor_pb_cpg <- cor(pubera_cpg, brevi_cpg, method = "spearman")
cor_pt_cpg <- cor(pubera_cpg, tenuis_cpg, method = "spearman")
cor_bt_cpg <- cor(brevi_cpg, tenuis_cpg, method = "spearman")

# CHG
cor_pb_chg <- cor(pubera_chg, brevi_chg, method = "spearman")
cor_pt_chg <- cor(pubera_chg, tenuis_chg, method = "spearman")
cor_bt_chg <- cor(brevi_chg, tenuis_chg, method = "spearman")

# CHG
cor_pb_chh <- cor(pubera_chh, brevi_chh, method = "spearman")
cor_pt_chh <- cor(pubera_chh, tenuis_chh, method = "spearman")
cor_bt_chh <- cor(brevi_chh, tenuis_chh, method = "spearman")

