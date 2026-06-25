SPdir="/netscratch/dep_mercier/grp_marques/marques/Rhync_tenuis_pangenome_project/Pangenome_subphasing"

for f in $SPdir/*_phase-results/*.ltr.insert.data
do
  cut -f3 $f | tail -n+2 >> ltr_insert_from_SubPhaser.txt
done
