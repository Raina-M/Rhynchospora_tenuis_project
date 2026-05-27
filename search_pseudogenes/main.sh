#!bin/bash

WD="/work/path/detect_pseudogenes"


# -------- 1. Extract functional proteins --------
# DNA fragments annotated in Helixer with RNA coverage rate > 0 are considered as functional genes.
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$7}' ${WD}/R_tenuis_REC_diploid_helixer.gtf > ${WD}/R_tenuis_REC_diploid_helixer.bed

# Compute RNA coverage for each helixer gene
bedtools coverage -a ${WD}/R_tenuis_REC_diploid_helixer.bed \
                  -b $WD/STAR_REC_dip_Aligned.sortedByCoord.out.filtered.bam > $WD/helixer_ann_with_RNA_covered.bed

awk '$8>0 {print $1"\t"$2"\t"$3}' $WD/helixer_ann_with_RNA_covered.bed > ${WD}/functional_genes.bed

# Extract the protein sequences of functional genes
bedtools intersect -a ${WD}/R_tenuis_REC_diploid_helixer.gtf -b ${WD}/functional_genes.bed -u > ${WD}/R_tenuis_REC_diploid_helixer_functional_genes.gtf

gffread -g ${WD}/R_tenuis_REC_diploid.fasta \
        -y ${WD}/proteins_of_functional_genes.fa \
        ${WD}/R_tenuis_REC_diploid_helixer_functional_genes.gtf

# -------- 2. Mask genome repeats and functional genes --------
# mask repeats based on earlGrey annotation
H1_ann="/path/Rhynchospora_tenuis_REC_h1.filteredRepeats.bed"
H2_ann="/path/Rhynchospora_tenuis_REC_h2.filteredRepeats.bed"

# concatenate the bed files that need to be masked
awk '$8>0 {print $1"\t"$2"\t"$3}' $WD/functional_genes.bed > ${WD}/masked_regions.bed
cut -f1-3 ${H1_ann} >> ${WD}/masked_regions.bed
cut -f1-3 ${H2_ann} >> ${WD}/masked_regions.bed

# hard mask
bedtools maskfasta -fi ${WD}/R_tenuis_REC_diploid.fasta \
                   -bed ${WD}/masked_regions.bed \
                   -fo ${WD}/R_tenuis_REC_diploid_HardMasked.fasta

# -------- 3. Align proteins to the masked genome --------
/netscratch/dep_mercier/grp_marques/bin/miniprot-0.18/miniprot \
  -Iut32 --gff ${WD}/R_tenuis_REC_diploid_HardMasked.fasta  \
  ${WD}/proteins_of_functional_genes.fa \
  > ${WD}/protein_map_to_maskedGenome.gff

# 
bash $WD/run_pipeline.sh ${WD}/protein_map_to_maskedGenome.gff ${WD}/output


# Pseudogenes Identification:
# | Type          | Origin             | Genomic Characteristics          | RNA-seq Signature                 |
# | ------------- | ------------------ | -------------------------------- | --------------------------------- |
# | Processed     | Retrotransposition | No introns, poly-A tail remnants | Usually silent or low expression  |
# | Non-processed | Gene Duplication   | Intron-exon structure preserved  | Often "dead" due to mutations     |
# | Unitary       | Mutation in situ   | No duplicate exists elsewhere    | Residual, non-functional mRNA     |