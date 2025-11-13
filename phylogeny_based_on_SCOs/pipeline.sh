WD="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/compare_genomes/phylogeny/breviuscula_as_outgroup/"

orthoSeqs="/netscratch/dep_mercier/grp_marques/marques/Rhync_tenuis_pangenome_project/GENESPACE/orthofinder/Results_Jun19/Single_Copy_Orthologue_Sequences"

cd $WD

### Align SCO seqs in the same group ###

# add suffix string to breviuscula,JGV16, JGV17 and REC because they can have the same protein IDs
string1="breviuscula"
string2="JGV-16"
string3="JGV-17"
string4="REC"

for fasta in $orthoSeqs/*.fa
do
  fname=`basename $fasta`
  #
  echo "Processing $fasta"
    
  awk -v s1="$string1" -v s2="$string2" -v s3="$string3" -v s4="$string4" '
    BEGIN {header_count = 0}
    /^>/ {
        header_count++
        if (header_count == 2) {
            $0 = $0 "_" s1
        }
        else if (header_count == 15 || header_count == 16) {
            $0 = $0 "_" s2
        }
        else if (header_count == 17 || header_count == 18) {
            $0 = $0 "_" s3
        }
        else if (header_count == 19 || header_count == 20) {
            $0 = $0 "_" s4
        }
    }
    {print}
    ' "$fasta" > $WD/$fname
  
  mafft --auto --thread 8 --maxiterate 10 $WD/$fname > $WD/MAFFT_alns/${fname%.fa}_aln.fasta
  rm $WD/$fname
done


### Concatenate SCO alignments ###
python3.7 concatenate_aln.py


### Run RAxML model to construct phylogeny and estimate divergence time ###
conda activate /netscratch/dep_mercier/grp_marques/bin/conda-envs/envs/raxml-ng

# build tree
raxml-ng --all --msa concatenated.fasta --model partitions.txt --thread 64




