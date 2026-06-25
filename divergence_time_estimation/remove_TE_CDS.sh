WD="/working/path"

cd $WD

# input=$1
input="$WD/dN_dS_ratio/CoGe_Kn_Ks/reference_hap1_vs_hap2_helixer.txt"

# only keep the CDS pairs with identity over 80
grep -v "^#" $input | awk -F"\t" '{split($4,a,"\\|\\|"); split($8,b,"\\|\\|"); if(a[9]>80 && b[9]>80) print $0}' > KnKs_over_80perc_identity.txt
# 21708 -> 19334

# ---------- remove the CDS pairs overlapping TEs ---------- #

# 1. Extract the coordinates of CDS pairs that have over 80% identity
cut -f4 KnKs_over_80perc_identity.txt | awk -F"[|][|]" '{print "Chr"$1"\t"$2"\t"$3"\t"NR}' > REC_hap1_KsKn_tmp.bed

cut -f8 KnKs_over_80perc_identity.txt | awk -F"[|][|]" '{print "Chr"$1"\t"$2"\t"$3"\t"NR}' > REC_hap2_KsKn_tmp.bed


# 2. remove the regions that overlap TE annotations
TEh1="$WD/TEs/REC_h1/Rhynchospora_tenuis_REC_h1_EarlGrey/Rhynchospora_tenuis_REC_h1_summaryFiles/Rhynchospora_tenuis_REC_h1.filteredRepeats.bed"
TEh2="$WD/TEs/REC_h2/Rhynchospora_tenuis_REC_h2_EarlGrey/Rhynchospora_tenuis_REC_h2_summaryFiles/Rhynchospora_tenuis_REC_h2.filteredRepeats.bed"

#  extract line numbers that not overlap with any TE annotations 
bedtools subtract -A -a REC_hap1_KsKn_tmp.bed -b $TEh1 | cut -f4 > REC_hap1_nonOverlapped_tmp.list
bedtools subtract -A -a REC_hap2_KsKn_tmp.bed -b $TEh2 | cut -f4 > REC_hap2_nonOverlapped_tmp.list

# only keep lines that both haps not overlapping TEs
cat REC_hap?_nonOverlapped_tmp.list | sort -n | uniq -c | awk '$1==2 {print $2}' > REC_nonOverlapped_tmp.list
sed -n -f <(sed 's/$/p/' REC_nonOverlapped_tmp.list) KnKs_over_80perc_identity.txt > KnKs_over_80perc_identity_remove_TEs.txt
# 19334 -> 11880

rm ./*tmp*

cut -f1 KnKs_over_80perc_identity_remove_TEs.txt > REC_Ks_over_80perc_identity_remove_TEs.list