#!bin/bash

WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/DNA_methylation"
seqdir="/biodata/dep_mercier/grp_marques/marques/other_Rhynchospora_seqdata/Rhync_tenuis_Methyl-seq"
REF="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/gene_screen/RNA_seq_mapping/R_tenuis_ref_diploid.fasta"

# process methy-seq data by bismark
BM_suite="/netscratch/dep_mercier/grp_marques/bin/marques-envs/bismark/bin"


# Prepare genome
# genome was already indexed in ChIP seq analyses, so use the bowtie index from ChIP
$BM_suite/bismark_genome_preparation --parallel 8 /netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/ChIP_seq/genome_index/


cd $WD
# Alignment
zcat ${seqdir}/*_R1_001.fastq.gz | gzip > $WD/methy_seq_R1.fastq.gz
zcat ${seqdir}/*_R2_001.fastq.gz | gzip > $WD/methy_seq_R2.fastq.gz

# no-mix no-discordant by default
${BM_suite}/bismark --parallel 32 --local \
  --prefix R_tenuis_dip \
  --genome /netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/ChIP_seq/genome_index/ \
  -o $WD \
  -1 $WD/methy_seq_R1.fastq.gz \
  -2 $WD/methy_seq_R2.fastq.gz -o $WD 2>&1 | tee -a $WD/bismark.log


# Deduplication
${BM_suite}/deduplicate_bismark -p --output_dir $WD --parallel 32 \
  --bam $WD/R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.bam 2>&1 | tee -a $WD/bismark_dedup.log


# Extract context-dependent (CpG/CHG/CHH) methylation
${BM_suite}/bismark_methylation_extractor -p \
  -o $WD --parallel 32 --gzip \
  --comprehensive --bedGraph \
  --split_by_chromosome \
  R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplicated.bam 2>&1 | tee -a $WD/bismark_extract.log


# visulize report
${BM_suite}/bismark2report --dir $WD \
  --alignment_report $WD/R_tenuis_dip.methy_seq_R1_bismark_bt2_PE_report.txt \
  --dedup_report $WD/R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplication_report.txt \
  --splitting_report $WD/R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplicated_splitting_report.txt \
  --mbias_report $WD/R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplicated.M-bias.txt




# Convert to bedGraph
${BM_suite}/bismark2bedGraph --dir $WD/bedgraph/ -o CpG.bedgraph CpG_context_R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplicated.txt.gz
${BM_suite}/bismark2bedGraph --CX --dir $WD/bedgraph/ -o CHG.bedgraph CHG_context_R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplicated.txt.gz
${BM_suite}/bismark2bedGraph --CX --dir $WD/bedgraph/ -o CHH.bedgraph CHH_context_R_tenuis_dip.methy_seq_R1_bismark_bt2_pe.deduplicated.txt.gz


# separate hap1 and hap2
cd $WD/bedgraph
for methyl in CpG CHG CHH
do
  zcat ${methyl}.bedgraph.gz | awk '$1=="Chr1_h2" || $1=="Chr2_h1"' | sort -k1,1 -k2,2n > ${methyl}_hap1.bedgraph
  zcat ${methyl}.bedgraph.gz | awk '$1=="Chr1_h1" || $1=="Chr2_h2"' | sort -k1,1 -k2,2n > ${methyl}_hap2.bedgraph  
done


# convert to bigwig
awk '$1=="Chr1_h2" || $1=="Chr2_h1" {print $1"\t"$2}' ${REF}.fai > hap1.genome
awk '$1=="Chr1_h1" || $1=="Chr2_h2" {print $1"\t"$2}' ${REF}.fai > hap2.genome

/netscratch/dep_mercier/grp_marques/bin/bigWig/bedGraphToBigWig CpG_hap1.bedgraph hap1.genome CpG_hap1.bw
/netscratch/dep_mercier/grp_marques/bin/bigWig/bedGraphToBigWig CpG_hap2.bedgraph hap2.genome CpG_hap2.bw
/netscratch/dep_mercier/grp_marques/bin/bigWig/bedGraphToBigWig CHG_hap1.bedgraph hap1.genome CHG_hap1.bw
/netscratch/dep_mercier/grp_marques/bin/bigWig/bedGraphToBigWig CHG_hap2.bedgraph hap2.genome CHG_hap2.bw
/netscratch/dep_mercier/grp_marques/bin/bigWig/bedGraphToBigWig CHH_hap1.bedgraph hap1.genome CHH_hap1.bw
/netscratch/dep_mercier/grp_marques/bin/bigWig/bedGraphToBigWig CHH_hap2.bedgraph hap2.genome CHH_hap2.bw
  
