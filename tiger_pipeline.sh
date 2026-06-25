#!/bin/bash

# Adapted from:
# tigerrun.sh runs TIGER pipeline on cluster
# 23 September 2015, by Beth Rowan
# Or follow the described pipeline in GitHub page:
# https://github.com/Imoteph/TIGER_Whole-Genome_Genotyping-by-Sequencing

# specify path to Java installation and java memory
export JAVA_HOME=/usr/bin/java
export _JAVA_OPTIONS=-Xmx1024M

# specify variables
TIGER=/netscratch/dep_mercier/grp_marques/bin/TIGER

WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/3_run_TIGER"
SNP_DIR="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/2_F1plants_SNP_calling"
outDir=$WD

# complete markers
MARKCOMP=${WD}/../1_marker_selection/converted_variant_merged.txt
# filtered markers
MARKCORR=${WD}/../1_marker_selection/Rtenuis_filtered_markers_merged.txt
chrsizes="/netscratch/dep_mercier/grp_marques/mzhang/HiC_maps/Rhynchospora_tenuis/references/Rtenuis.hap1.hifiasm.phased.hic.chrom.sizes"

### ---------- Step 1: Create input files ----------
# marker files need to be modified before being input to get_subset.pl
# because column names of input.out and marker files are different
cd $outDir
sed 's/Chr1_h1/1/g
     s/Chr2_h1/2/g' < $MARKCOMP > complete_markers.txt
sed 's/Chr1_h1/1/g
     s/Chr2_h1/2/g' < $MARKCORR > corrected_markers.txt

# Input files can be generated from either sam (samtotiger.sh) or bam (bamtotiger.sh) files.
# but we already have snp files (VCF), so only adopt part of the script to generate files we need
cd $SNP_DIR
for snpfile in *.vcf
do
	# create a single directory for each sample, named by sample names
	mkdir ${outDir}/${snpfile%_*}
	cd ${outDir}/${snpfile%_*}

	#extract DP4 lines
	bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ${SNP_DIR}/$snpfile -o ${snpfile%_*}.comma.txt

	#replace comma separators with tabs
	tr ',' '\t' < ${snpfile%_*}.comma.txt > ${snpfile%_*}.tabbed.txt
	rm ${snpfile%_*}.comma.txt

	# get rid of mito and cp reads;
	# create columns for read counts for ref allele and alt allele by adding the first two of the DP4 fields and the second two of the DP4 fields.
	# INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
	awk '{if ($1 !="ChrC" && $1!="ChrM") print $1 "\t" $2 "\t" $3 "\t" $8+$9 "\t"  $4 "\t"  $10+$11}' \
	    ${snpfile%_*}.tabbed.txt > ${snpfile%_*}.input.temp
	rm ${snpfile%_*}.tabbed.txt

	# change chromosome column to numbers only (to match TIGER input)
	sed 's/Chr1_h1/1/g
             s/Chr2_h1/2/g' < ${snpfile%_*}.input.temp > input.out
	#bash $TIGER/chr_convert.sh ${snpfile%_*}.input.temp
	rm ${snpfile%_*}.input.temp

	# filter indels only (to generate "complete" file for TIGER)
	perl $TIGER/get_subset.pl input.out 1,2 ${outDir}/complete_markers.txt 1,2 0 > input_complete.txt

	# use stricter filters (to generate "corrected" file for TIGER)
	perl $TIGER/get_subset.pl input.out 1,2 ${outDir}/corrected_markers.txt 1,2 0 > input_corrected.txt

	rm input.out


### ---------- Step 2: Run base caller on corrected input file ----------
	# -n specifies biparental
	java -jar $TIGER/base_caller.jar -r input_corrected.txt -o allele_count_base_call.txt -n bi
	# output file is a single line for each chromosome with all of the positions assigned a base caller genotype
	# genotype notation uses C and L for the two different parental alleles


### ---------- Step 3: Run allele frequency estimator ----------
        # apply a sliding window to estimates the local allele frequency for N (-w argument) adjacent markers
	java -jar $TIGER/allele_freq_estimator.jar -r input_corrected.txt -o frequencies_for_bmm.txt -n bi -w 1000


### ---------- Step 4: Apply the beta mixture model ----------
	R --slave --vanilla --args frequencies_for_bmm.txt bmm.intersections.txt < $TIGER/beta_mixture_model.R
	#Or: Rscript --vanilla $TIGER/beta_mixture_model.R frequencies_for_bmm.txt bmm.intersections.txt


### ---------- Step 5: Prepare files for HMM probability estimation ---------- 
	# using the output from BASECALLER and beta mixiture model
	perl $TIGER/prep_prob.pl -s ${snpfile%_*} -m input_corrected.txt -b allele_count_base_call.txt -c $chrsizes -o file_for_probabilities.txt


### ---------- Step 6: Calculate transmission and emission probabilities for the HMM General command ----------
	# This gives two output files: 1. sample_hmm_model (probabilities for the HMM)
	#			       2. sample_sliding_window (genotyping only based on the sliding window)
	perl $TIGER/hmm_prob.pl -s frequencies_for_bmm.txt -p file_for_probabilities.txt -o ${snpfile%_*} -a bmm.intersections.txt -c $chrsizes


### ---------- Step 7: Run the HMM ----------

	java -jar $TIGER/hmm_play.jar -r allele_count_base_call.txt -o hmm.out.txt -t bi -z ${snpfile%_*}_hmm_model


### ---------- Step 8: Get rough estimate of recombination breakpoint positions ----------
	perl $TIGER/prepare_break.pl -s ${snpfile%_*} -m input_corrected.txt -b hmm.out.txt -c $chrsizes -o rough_COs.txt


### ---------- Step 9: Refine recombination breaks ---------- 
	# using nearby markers that had been previously filtered out
	perl $TIGER/refine_recombination_break.pl input_complete.txt rough_COs.breaks.txt


### ---------- Step 10: Smooth out breaks ----------
	perl $TIGER/breaks_smoother.pl -b rough_COs.refined.breaks.txt -o ${snpfile%_*}.corrected.refined.breaks.txt


### ---------- Visualize output ----------
	sh ${WD}/visualize_result.sh ${snpfile%_*} $outDir $chrsizes

done
