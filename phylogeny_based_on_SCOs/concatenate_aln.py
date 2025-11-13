#!/usr/bin/env python3
from Bio import SeqIO
import glob

# === USER SETTINGS ===
input_pattern = "MAFFT_alns/*_aln.fasta"   # path to your MAFFT outputs
output_fasta = "concatenated.fasta"
output_part = "partitions.txt"

# canonical species order & names
species_order = [
    "RHYNC_austrobrasiliensis",
    "RHYNC_breviuscula",
    "RHYNC_tenuis_PECP-36_h1",
    "RHYNC_tenuis_PECP-36_h2",
    "RHYNC_tenuis_JGV-89_h1",
    "RHYNC_tenuis_JGV-89_h2",
    "RHYNC_tenuis_PECP-47_h1",
    "RHYNC_tenuis_PECP-47_h2",
    "RHYNC_tenuis_PECP-48_h1",
    "RHYNC_tenuis_PECP-48_h2",
    "RHYNC_tenuis_PECP-35-6_h1",
    "RHYNC_tenuis_PECP-35-6_h2",
    "RHYNC_tenuis_PECP-36-7_h1",
    "RHYNC_tenuis_PECP-36-7_h2",
    "RHYNC_tenuis_JGV-16_h1",
    "RHYNC_tenuis_JGV-16_h2",
    "RHYNC_tenuis_JGV-17_h1",
    "RHYNC_tenuis_JGV-17_h2",
    "RHYNC_tenuis_REC_h1",
    "RHYNC_tenuis_REC_h2"
]

# === SCRIPT BEGINS ===
files = sorted(glob.glob(input_pattern))
if not files:
    raise SystemExit("No alignment files found — check your input path.")

print(f"Found {len(files)} alignment files.")
n_species = len(species_order)

# Initialize concatenated sequences
concat = {sp: "" for sp in species_order}
partitions = []
pos = 1

for aln_file in files:
    records = list(SeqIO.parse(aln_file, "fasta"))
    if len(records) != n_species:
        raise ValueError(
            f"{aln_file} has {len(records)} sequences but expected {n_species}."
        )

    aln_len = len(records[0].seq)
    gene = aln_file.split("/")[-1].split(".")[0]

    print(f"Adding {gene}: length {aln_len}")
    for i, sp in enumerate(species_order):
        seq = str(records[i].seq).replace("*", "-").upper()
        concat[sp] += seq

    partitions.append(f"LG+G+F, {gene} = {pos}-{pos + aln_len - 1}")
    pos += aln_len

# Write concatenated FASTA
with open(output_fasta, "w") as out:
    for sp in species_order:
        out.write(f">{sp}\n{concat[sp]}\n")

# Write partitions
with open(output_part, "w") as part:
    part.write("\n".join(partitions))

print(f"\nWrote {output_fasta} and {output_part}")
print("Each sequence length:", pos - 1)

