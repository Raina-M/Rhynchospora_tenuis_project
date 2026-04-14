
WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/DNA_methylation/metaplots"

module load deeptools/v3.5.5

cd $WD

# REC haplotype 1 Tyba arrays
grep "TRC_1;" /netscratch/dep_mercier/grp_marques/Laura/Rhynchospora/Pangenome_Rtenuis/Rhync_tenuis_REC-hap1/TideCluster/Rhync_tenuis_REC-hap1__TideCluster/tidecluster_annotation.gff3 | cut -f1,4,5 > REC_hap1_Tyba_arrays.bed

grep "TRC_1;" /netscratch/dep_mercier/grp_marques/Laura/Rhynchospora/Pangenome_Rtenuis/Rhync_tenuis_REC-hap2/TideCluster/Rhync_tenuis_REC-hap2__TideCluster/tidecluster_annotation.gff3 | cut -f1,4,5 > REC_hap2_Tyba_arrays.bed


# REC hap1 genes
H1_genes="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_REC.hap1.chr.helixer.gff3"
# convert to bed format
grep -v "^#" $H1_genes | awk '$3=="gene" {print $1"\t"$4"\t"$5"\t"$7}' > $WD/H1_genes.bed

H2_genes="/biodata/dep_mercier/grp_marques/marques/Assemblies/Rhynch_tenuis_pangenome_assemblies/Rhync_tenuis_REC.hap2.chr.helixer.gff3"
# convert to bed format
grep -v "^#" $H2_genes | awk '$3=="gene" {print $1"\t"$4"\t"$5"\t"$7}' > $WD/H2_genes.bed



# TEs
H1_TE="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/compare_genomes/TEs/REC_h1/Rhynchospora_tenuis_REC_h1_EarlGrey/Rhynchospora_tenuis_REC_h1_summaryFiles/Rhynchospora_tenuis_REC_h1.filteredRepeats.bed"
H2_TE="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/compare_genomes/TEs/REC_h2/Rhynchospora_tenuis_REC_h2_EarlGrey/Rhynchospora_tenuis_REC_h2_summaryFiles/Rhynchospora_tenuis_REC_h2.filteredRepeats.bed"


# remove genes that overlap with TEs
bedtools subtract -A -a $WD/H1_genes.bed -b $H1_TE > $WD/H1_genes_purged.bed
bedtools subtract -A -a $WD/H2_genes.bed -b $H2_TE > $WD/H2_genes_purged.bed


############# CENH3 Metaplot over TYBA ############
# CpG
computeMatrix scale-regions \
  -S $WD/../bedgraph/CpG_hap1.bw $WD/../bedgraph/CHG_hap1.bw $WD/../bedgraph/CHH_hap1.bw\
  --regionsFileName REC_hap1_Tyba_arrays.bed H1_genes_purged.bed $H1_TE \
  -a 2000 -b 2000 --numberOfProcessors 40 \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 4000 \
  --afterRegionStartLength 2000 \
  --smartLabels \
  -o REC_H1_DNAmethyl_matrix.gz

plotHeatmap -m REC_H1_DNAmethyl_matrix.gz \
            --plotFileFormat "pdf" -o REC_H1_DNAmethyl_metaplot_rescale.pdf \
            -y log2RPKM \
            --startLabel "Start" --endLabel "End" \
            --heatmapWidth 12



# haplotype 2
computeMatrix scale-regions \
  -S $WD/../bedgraph/CpG_hap2.bw $WD/../bedgraph/CHG_hap2.bw $WD/../bedgraph/CHH_hap2.bw\
  --regionsFileName REC_hap2_Tyba_arrays.bed H2_genes_purged.bed $H2_TE \
  -a 2000 -b 2000 --numberOfProcessors 40 \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 4000 \
  --afterRegionStartLength 2000 \
  --smartLabels \
  -o REC_H2_DNAmethyl_matrix.gz
  

plotHeatmap -m REC_H2_DNAmethyl_matrix.gz \
            --plotFileFormat "pdf" -o REC_H2_DNAmethyl_metaplot.pdf \
            -y log2RPKM \
            --heatmapWidth 12 \
            --startLabel "Start" --endLabel "End"

