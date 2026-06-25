#!/bin/bash
#SBATCH -p normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-type=NONE

ref=$1
barcode=$2
markers=$3

TIGER=/netscratch/dep_mercier/grp_marques/bin/TIGER


##############
# SNP calling
##############
echo "mpileup to generating genotype likelihood"
bcftools mpileup -Oz --threads 1 --min-MQ 1 -B --annotate AD -o ${barcode}.mpileup.gz -f $ref ${barcode}.sorted.bam

echo -e "*\t*\t*\t*\t2" > ploidy.txt

echo "call variants"
bcftools call -Am -Oz --threads 1 --ploidy-file ploidy.txt \
	      ${barcode}.mpileup.gz -o ${barcode}.vcf.gz

####################################################
# Convert variants to desired format for CO calling
####################################################
# extract DP4 lines
bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ${barcode}.vcf.gz -o comma.txt

# replace comma separators with tabs
tr ',' '\t' < comma.txt > tabbed.txt
rm comma.txt

awk '{print $1 "\t" $2 "\t" $3 "\t" $8+$9 "\t"  $4 "\t"  $10+$11}' \
    tabbed.txt > input.tmp
rm tabbed.txt

# change chromosome column to numbers only (to match TIGER input)
sed 's/Chr1_h1/1/g
     s/Chr2_h1/2/g' < input.tmp > input.out

rm input.tmp

# use stricter filters (to generate "corrected" file for TIGER)
perl $TIGER/get_subset.pl input.out 1,2 $markers 2,3 0 > input_corrected.txt

rm input.out
