import gzip
import numpy as np

# Input file from command line
file_path="/netscratch/dep_mercier/grp_marques/mzhang/Rtenuis_CO_calling/DNA_methylation/metaplots/REC_H1_DNAmethyl_matrix.gz"

# Matrix header of R. pubera DNA methylation data
#@{"upstream":[2000,2000,2000],
#  "downstream":[2000,2000,2000],
#  "body":[4000,4000,4000],
#  "bin size":[10,10,10],
#  "ref point":[null,null,null],
#  "verbose":false,
#  "bin avg type":"mean",
#  "missing data as zero":false,
#  "min threshold":null,
#  "max threshold":null,
#  "scale":1,
#  "skip zeros":false,
#  "nan after end":false,
#  "proc number":40,
#  "sort regions":"keep",
#  "sort using":"mean",
#  "unscaled 5 prime":[0,0,0],
#  "unscaled 3 prime":[0,0,0],
#  "group_labels":["REC_hap1_Tyba_arrays.bed","H1_genes_purged.bed","Rhynchospora_tenuis_REC_h1.filteredRepeats.bed"],
#  "group_boundaries":[0,1282,18729,205319],
#  "sample_labels":["CpG_hap1","CHG_hap1","CHH_hap1"],
#  "sample_boundaries":[0,800,1600,2400]}


# Based on on matrix header, here we calculate the column means for the Tyba arrays only
START_LINE = 0
END_LINE = 1282

# The boundaries for the 4 ChIP signals (each has 800 bins)
samples = {
    "CpG":    (0, 800),
    "CHG":  (800, 1600),
    "CHH":  (1600, 2400)
}

# Arrays to hold the sum and counts for the 2400 columns
sums = np.zeros(2400)
counts = np.zeros(2400)

print(f"Streaming {file_path}...")
print(f"Extracting lines {START_LINE} to {END_LINE} (Tyba Arrays)...")

with gzip.open(file_path, 'rt') as f:
    for line_num, line in enumerate(f, start=1):
        
        # Stop reading entirely once we pass the Tyba section (Saves massive time)
        if line_num > END_LINE:
            break
            
        # Only process lines within the Tyba boundaries
        if line_num >= START_LINE:
            # Split the line, skip the first 6 metadata columns
            vals = line.split('\t')[6:]
            row = np.array(vals, dtype=float)
            
            # Find non-nan values and add them to our running totals
            mask = ~np.isnan(row)
            sums[mask] += row[mask]
            counts[mask] += 1

print("Calculating means...")
# Calculate the final means (ignoring division by zero if a column is 100% nans)
with np.errstate(invalid='ignore'):
    means = np.divide(sums, counts, out=np.zeros_like(sums), where=counts!=0)

# Split the 3200 means into 4 files
for name, (start, end) in samples.items():
    sample_means = means[start:end]
    out_file = f"{name}_Tyba_colMeans.txt"
    
    with open(out_file, "w") as out:
        out.write("\t".join(map(str, sample_means)) + "\n")
        
    print(f"Saved {name} means ({len(sample_means)} bins) to -> {out_file}")

print("Done!")
