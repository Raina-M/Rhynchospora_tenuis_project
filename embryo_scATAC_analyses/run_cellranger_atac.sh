#!bin/bash

WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/BD_scRNA/ref_embryo_10XATAC"

REF="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/scATACseq/cellranger_atac/references/rtenuis_hap1_chrs"

seqDir="/biodata/dep_mercier/grp_marques/marques/Hologen/SINGLE_CELL_SEQ_DATA/7248_scATAC_seq/"

markers="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/sc_COs_manual/corrected_markers.txt"

##### cellranger-atac #####
cd $WD
cellranger-atac count --id=R_tenuis \
                      --reference=$REF \
                      --fastqs=$seqDir  \
                      --sample=7248B  \
                      --localcores=32 \
                      --localmem=64


##### 2. split bam from cellranger output #####
BAM="${WD}/R_tenuis/outs/possorted_bam.bam"
samtools view  $BAM | cut -f 12- | tr "\t" "\n"  | grep  "^CB:Z:"  | cut -d ':' -f 3 | sort --parallel 16 | uniq -c > ${WD}/read_num.stats
# 540 cells with >=10000 reads, account for 739341824/804626684 (92%) of all reads


##############
# Demultiplex
##############
# exclude unmapped reads to reduce running time
samtools view -hb -F4 -@16 $BAM > $WD/aln_tmp.bam

awk '$1>=10000 {print $2}' $WD/read_num.stats | while read barcode
do  
  [ ! -d ${WD}/demultiplex/${barcode} ] && mkdir -p "${WD}/demultiplex/$barcode"
  cd ${WD}/demultiplex/$barcode
  
  sbatch $WD/scripts/slurm_demux.sh $WD/aln_tmp.bam $barcode $WD/read_num.stats $WD
done

rm $WD/aln_tmp.bam

# due to the low number of cells, no filtering at this stage


#####################################
# SNP calling and convert to markers
#####################################
awk '$1>=10000 {print $2}' $WD/read_num.stats | while read barcode
do
  cd ${WD}/demultiplex/$barcode
  sbatch $WD/scripts/slurm_snpCall_and_conversion.sh $REF/fasta/genome.fa $barcode $markers
done


##############################
# Count markers and filtering
##############################
awk '$1>=10000 {print $2}' $WD/read_num.stats | while read barcode
do
	cd $WD/demultiplex/$barcode
    awk '$4!=0 || $6!=0' input_corrected.txt > input_corrected_tmp.txt
    marker_num=`wc -l input_corrected_tmp.txt | cut -d" " -f1`
    # count genotype switching times
    awk '{
        if      ( $4 / ($4+$6) < 0.2) print 0;
        else if ( $4 / ($4+$6) > 0.8) print 1;
        else print 0.5;
    }' input_corrected_tmp.txt > smt_genotypes_tmp.txt
    head -n-1 smt_genotypes_tmp.txt > tmp1.txt
    tail -n+2 smt_genotypes_tmp.txt > tmp2.txt
    switch_num=`paste tmp1.txt tmp2.txt | awk '$2-$1!=0' | wc -l`
    echo -e $barcode"\t"$marker_num"\t"$switch_num >> $WD/switches.stats 

    rm *tmp*
done


#############
# Genotyping
#############
awk '$2>500 {print $1}' $WD/switches.stats | while read barcode
do
  Rscript $WD/scripts/hapCO_identification.R -i $WD/demultiplex/$barcode/input_corrected.txt -p ${barcode} -g $REF/fasta/genome.fa.fai -o ${WD}/CO_results_3M10 -c 500 -s 3000000 -n 10
done


### nomalized genotyping
# Create sliding windows
awk '{print $1"\t"$2}' $REF/fasta/genome.fa.fai | sed 's/Chr1_h1/1/g;s/Chr2_h1/2/g' > $WD/R_tenuis.genome
bedtools makewindows -g $WD/R_tenuis.genome -w 1000000 -s 500000 > $WD/windows_1M_500k.bed 


# Calculate normAF per window and then genotyping
awk '$2>500 {print $1}' $WD/switches.stats | while read barcode
do
  # normalized AF
  awk '$4!=0 || $6!=0 {print $1"\t"$2"\t"$2+1"\t"2*$4/($4+$6)-1}' $WD/demultiplex/$barcode/input_corrected.txt > $WD/${barcode}_normAF.bed
  
  # calculate normalized AF in windows
  bedtools map -a $WD/windows_1M_500k.bed -b $WD/${barcode}_normAF.bed -c 4 -o mean -null 0 > $WD/${barcode}_normAF_by_win.bed
  
  echo "Start to genotype $barcode ..."
  Rscript $WD/scripts/genotype_conversion.R $WD/${barcode}_normAF_by_win.bed $REF/fasta/genome.fa.fai ${WD}/genotype_norm/${barcode}_norm_genotype.pdf
  
  # clean up
  rm $WD/${barcode}_normAF.bed
  rm $WD/${barcode}_normAF_by_win.bed
done
