#!/bin/bash
#SBATCH -p normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --mail-type=NONE

input=$1
barcode=$2
readNum=$3
outDir=$4

# extract read names with this barcode
samtools view -h $input | egrep "^@|CB:Z:${barcode}" | samtools sort -@4 -o ${barcode}.sorted.bam

# total read number
read_num=`awk -v bc=$barcode '$2==bc {print $1}' $readNum`
# unique mapping read number
uniq_map=`samtools view ${barcode}.sorted.bam | awk '$5==60' | wc -l`

echo $barcode" "$read_num" "$uniq_map >> $outDir/aligned_num.stats


