#!/bin/sh
WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/compare_genomes/SYRI"

PY=/netscratch/dep_mercier/grp_marques/bin/marques-envs/syri/bin/python3.8
PLOTsr=/netscratch/dep_mercier/grp_marques/bin/marques-envs/syri/bin/plotsr

# genomes
REF1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_ref.hap1.chr.fasta"
REF2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_ref.hap2.chr.fasta"
JGV016_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6523A.hap1.chr.fasta"
JGV016_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6523A.hap2.chr.fasta"
JGV017_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6524A.hap1.chr.fasta"
JGV017_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6524A.hap2.chr.fasta"
JGV_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344A.JGV.hap1.chr.fasta"
JGV_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344A.JGV.hap2.chr.fasta"
PECP2_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344B.PECP2.hap1.chr.fasta"
PECP2_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344B.PECP2.hap2.chr.fasta"
PECP_035_6_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6365A.035_5.hap1.chr.fasta"
PECP_035_6_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6365A.035_5.hap2.chr.fasta"
PECP3_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344C.PECP3.hap1.chr.fasta"
PECP3_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6344C.PECP3.hap2.chr.fasta"
PECP_036_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6228D.036.hap1.chr.fasta"
PECP_036_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6228D.036.hap2.chr.fasta"
PECP_036_7_1="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6365B.036_7.hap1.chr.fasta"
PECP_036_7_2="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_6365B.036_7.hap2.chr.fasta"



######### mapping #########
cd $WD

# REF_h1 vs. REF_h2
minimap2 -ax asm5 --eqx -t 16 $REF1 $REF2 | samtools sort -@8 > $WD/REF1_REF2.sorted.bam
samtools index -@4 $WD/REF1_REF2.sorted.bam

# REF_h2 vs. JGV016_h1
minimap2 -ax asm5 --eqx -t 16 $REF2 $JGV016_1 | samtools sort -@8 > $WD/REF-h2_JGV016-h1.sorted.bam
samtools index -@4 $WD/REF-h2_JGV016-h1.sorted.bam

# JGV016_h1 vs. JGV016_h2
minimap2 -ax asm5 --eqx -t 16 $JGV016_1 $JGV016_2 | samtools sort -@8 > $WD/JGV016-h1_JGV016-h2.sorted.bam
samtools index -@4 $WD/JGV016-h1_JGV016-h2.sorted.bam

# JGV016_h2 vs. JGV017_h1
minimap2 -ax asm5 --eqx -t 16 $JGV016_2 $JGV017_1 | samtools sort -@8 > $WD/JGV016-h2_JGV017-h1.sorted.bam
samtools index -@4 $WD/JGV016-h2_JGV017-h1.sorted.bam

# JGV017_h1 vs. JGV017_h2
minimap2 -ax asm5 --eqx -t 16 $JGV017_1 $JGV017_2 | samtools sort -@8 > $WD/JGV017-h1_JGV017-h2.sorted.bam
samtools index -@4 $WD/JGV017-h1_JGV017-h2.sorted.bam

# JGV017_h2 vs. JGV_h1
minimap2 -ax asm5 --eqx -t 16 $JGV017_2 $JGV_1 | samtools sort -@8 > $WD/JGV017-h2_JGV-h1.sorted.bam
samtools index -@4 $WD/JGV017-h2_JGV-h1.sorted.bam

# JGV_h1 vs. JGV_h2
minimap2 -ax asm5 --eqx -t 16 $JGV_1 $JGV_2 | samtools sort -@8 > $WD/JGV-h1_JGV-h2.sorted.bam
samtools index -@4 $WD/JGV-h1_JGV-h2.sorted.bam

# JGV_h2 vs. PECP2_h1
minimap2 -ax asm5 --eqx -t 16 $JGV_2 $PECP2_1 | samtools sort -@8 > $WD/JGV-h2_PECP2-h1.sorted.bam
samtools index -@4 $WD/JGV-h2_PECP2-h1.sorted.bam

# PECP2_h1 vs. PECP2_h2
minimap2 -ax asm5 --eqx -t 16 $PECP2_1 $PECP2_2 | samtools sort -@8 > $WD/PECP2-h1_PECP2-h2.sorted.bam
samtools index -@4 $WD/PECP2-h1_PECP2-h2.sorted.bam

# PECP2_h2 vs. 035-6_h1
minimap2 -ax asm5 --eqx -t 16 $PECP2_2 $PECP_035_6_1 | samtools sort -@8 > $WD/PECP2-h2_035-6-h1.sorted.bam
samtools index -@4 $WD/PECP2-h2_035-6-h1.sorted.bam

# 035-6_h1 vs. 035-6_h2
minimap2 -ax asm5 --eqx -t 16 $PECP_035_6_1 $PECP_035_6_2 | samtools sort -@8 > $WD/035-6-h1_035-6-h2.sorted.bam
samtools index -@4 $WD/035-6-h1_035-6-h2.sorted.bam

# 035-6_h2 vs. PECP3_h1
minimap2 -ax asm5 --eqx -t 16 $PECP_035_6_2 $PECP3_1 | samtools sort -@8 > $WD/035-6-h2_PECP3-h1.sorted.bam
samtools index -@4 $WD/035-6-h2_PECP3-h1.sorted.bam

# PECP3_h1 vs. PECP3_h2
minimap2 -ax asm5 --eqx -t 16 $PECP3_1 $PECP3_2 | samtools sort -@8 > $WD/PECP3-h1_PECP3-h2.sorted.bam
samtools index -@4 $WD/PECP3-h1_PECP3-h2.sorted.bam

# PECP3_h2 vs. 036_h1
minimap2 -ax asm5 --eqx -t 16 $PECP3_2 $PECP_036_1 | samtools sort -@8 > $WD/PECP3-h2_036-h1.sorted.bam
samtools index -@4 $WD/PECP3-h2_036-h1.sorted.bam

# 036_h1 vs. 036_h2
minimap2 -ax asm5 --eqx -t 16 $PECP_036_1 $PECP_036_2 | samtools sort -@8 > $WD/036-h1_036-h2.sorted.bam
samtools index -@4 $WD/036-h1_036-h2.sorted.bam

# 036_h2 vs. 036-7_h1
minimap2 -ax asm5 --eqx -t 16 $PECP_036_2 $PECP_036_7_1 | samtools sort -@8 > $WD/036-h2_036-7-h1.sorted.bam
samtools index -@4 $WD/036-h2_036-7-h1.sorted.bam

# 036-7_h1 vs. 036-7_h2
minimap2 -ax asm5 --eqx -t 16 $PECP_036_7_1 $PECP_036_7_2 | samtools sort -@8 > $WD/036-7-h1_036-7-h2.sorted.bam
samtools index -@4 $WD/036-7-h1_036-7-h2.sorted.bam



######### syri #########
conda activate /netscratch/dep_mercier/grp_marques/bin/marques-envs/syri

# --- REF1 vs. REF2 ---
mkdir $WD/1_REF1_REF2
cd $WD/1_REF1_REF2
# run syri
syri -c $WD/REF1_REF2.sorted.bam -r $REF1 -q $REF2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- REF2 vs. JGV016-h1 ---
mkdir $WD/2_REF-h2_JGV016-h1
cd $WD/2_REF-h2_JGV016-h1
# run syri
syri -c $WD/REF-h2_JGV016-h1.sorted.bam -r $REF2 -q $JGV016_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- JGV016-h1 vs. JGV016-h2 ---
mkdir $WD/3_JGV016-h1_JGV016-h2
cd $WD/3_JGV016-h1_JGV016-h2
# run syri
syri -c $WD/JGV016-h1_JGV016-h2.sorted.bam -r $JGV016_1 -q $JGV016_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- JGV016-h2 vs. JGV017-h1 ---
mkdir $WD/4_JGV016-h2_JGV017-h1
cd $WD/4_JGV016-h2_JGV017-h1
# run syri
syri -c $WD/JGV016-h2_JGV017-h1.sorted.bam -r $JGV016_2 -q $JGV017_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- JGV017-h1 vs. JGV017-h2 ---
mkdir $WD/5_JGV017-h1_JGV017-h2
cd $WD/5_JGV017-h1_JGV017-h2
# run syri
syri -c $WD/JGV017-h1_JGV017-h2.sorted.bam -r $JGV017_1 -q $JGV017_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- JGV017-h2 vs. JGV-h1 ---
mkdir $WD/6_JGV017-h2_JGV-h1
cd $WD/6_JGV017-h2_JGV-h1
# run syri
syri -c $WD/JGV017-h2_JGV-h1.sorted.bam -r $JGV017_2 -q $JGV_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- JGV-h1 vs. JGV-h2 ---
mkdir $WD/7_JGV-h1_JGV-h2
cd $WD/7_JGV-h1_JGV-h2
# run syri
syri -c $WD/JGV-h1_JGV-h2.sorted.bam -r $JGV_1 -q $JGV_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- JGV-h2 vs. PECP2-h1 ---
mkdir $WD/8_JGV-h2_PECP2-h1
cd $WD/8_JGV-h2_PECP2-h1
# run syri
syri -c $WD/JGV-h2_PECP2-h1.sorted.bam -r $JGV_2 -q $PECP2_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- PECP2-h1 vs. PECP2-h2 ---
mkdir $WD/9_PECP2-h1_PECP2-h2
cd $WD/9_PECP2-h1_PECP2-h2
# run syri
syri -c $WD/PECP2-h1_PECP2-h2.sorted.bam -r $PECP2_1 -q $PECP2_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- PECP2-h2 vs. 035-6-h1 ---
mkdir $WD/10_PECP2-h2_035-6-h1
cd $WD/10_PECP2-h2_035-6-h1
# run syri
syri -c $WD/PECP2-h2_035-6-h1.sorted.bam -r $PECP2_2 -q $PECP_035_6_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- 035-6-h1 vs. 035-6-h2 ---
mkdir $WD/11_035-6-h1_035-6-h2
cd $WD/11_035-6-h1_035-6-h2
# run syri
syri -c $WD/035-6-h1_035-6-h2.sorted.bam -r $PECP_035_6_1 -q $PECP_035_6_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- 035-6-h2 vs. PECP3-h1 ---
mkdir $WD/12_035-6-h2_PECP3-h1
cd $WD/12_035-6-h2_PECP3-h1
# run syri
syri -c $WD/035-6-h2_PECP3-h1.sorted.bam -r $PECP_035_6_2 -q $PECP3_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- PECP3-h1 vs. PECP3-h2 ---
mkdir $WD/13_PECP3-h1_PECP3-h2
cd $WD/13_PECP3-h1_PECP3-h2
# run syri
syri -c $WD/PECP3-h1_PECP3-h2.sorted.bam -r $PECP3_1 -q $PECP3_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- PECP3-h2 vs. 036-h1 ---
mkdir $WD/14_PECP3-h2_036-h1
cd $WD/14_PECP3-h2_036-h1
# run syri
syri -c $WD/PECP3-h2_036-h1.sorted.bam -r $PECP3_2 -q $PECP_036_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- 036-h1 vs. 036-h2 ---
mkdir $WD/15_036-h1_036-h2
cd $WD/15_036-h1_036-h2
# run syri
syri -c $WD/036-h1_036-h2.sorted.bam -r $PECP_036_1 -q $PECP_036_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- 036-h2 vs. 036-7-h1 ---
mkdir $WD/16_036-h2_036-7-h1
cd $WD/16_036-h2_036-7-h1
# run syri
syri -c $WD/036-h2_036-7-h1.sorted.bam -r $PECP_036_2 -q $PECP_036_7_1 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

# --- 036-7-h1 vs. 036-7-h2 ---
mkdir $WD/17_036-7-h1_036-7-h2
cd $WD/17_036-7-h1_036-7-h2
# run syri
syri -c $WD/036-7-h1_036-7-h2.sorted.bam -r $PECP_036_7_1 -q $PECP_036_7_2 -k -F B --nc 16
# plot
$PY $PLOTsr --sr syri.out --genomes genomes.txt -o syri_50k.pdf -S 0.7 -W 6 -H 4 -f 10 --itx -s 50000

######### plot #########
$PY $PLOTsr --sr $WD/1_REF1_REF2/syri.out \
            --sr $WD/2_REF-h2_JGV016-h1/syri.out \
            --sr $WD/3_JGV016-h1_JGV016-h2/syri.out \
            --sr $WD/4_JGV016-h2_JGV017-h1/syri.out \
            --sr $WD/5_JGV017-h1_JGV017-h2/syri.out \
            --sr $WD/6_JGV017-h2_JGV-h1/syri.out \
            --sr $WD/7_JGV-h1_JGV-h2/syri.out \
            --sr $WD/8_JGV-h2_PECP2-h1/syri.out \
            --sr $WD/9_PECP2-h1_PECP2-h2/syri.out \
            --sr $WD/10_PECP2-h2_035-6-h1/syri.out \
            --sr $WD/11_035-6-h1_035-6-h2/syri.out \
            --sr $WD/12_035-6-h2_PECP3-h1/syri.out \
            --sr $WD/13_PECP3-h1_PECP3-h2/syri.out \
            --sr $WD/14_PECP3-h2_036-h1/syri.out \
            --sr $WD/15_036-h1_036-h2/syri.out \
            --sr $WD/16_036-h2_036-7-h1/syri.out \
            --sr $WD/17_036-7-h1_036-7-h2/syri.out \
            --genomes $WD/all_genomes.txt \
            -o $WD/syri_all_haps_50k_with_markers_newColors.pdf \
            -S 0.7 -W 5 -H 15 -f 10 \
            --itx -s 50000 \
            --markers 35S_rDNA_markers.txt

