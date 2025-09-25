import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# Parameters
n_permutations = 100  # Increase permutation number for smoother distributions

# Load Orthogroups.GeneCount.tsv
df = pd.read_csv("Orthogroups_GeneCount_all_tenius_haps.tsv", sep='\t', index_col=0)

# List of genome names (individuals)
genomes = list(df.columns)
n_genomes = len(genomes)

# Storage for all permutations
pan_all = np.zeros((n_permutations, n_genomes))
core_all = np.zeros((n_permutations, n_genomes))

for p in range(n_permutations):
    random.shuffle(genomes)
    pan_tmp = []
    core_tmp = []

    for i in range(n_genomes):
        subset = df[genomes[:i+1]]
        present_in_subset = subset[(subset > 0).any(axis=1)].index
        present_in_all = subset[(subset > 0).all(axis=1)].index

        pan_tmp.append(len(present_in_subset))
        core_tmp.append(len(present_in_all))

    pan_all[p, :] = pan_tmp
    core_all[p, :] = core_tmp

# Create boxplots for each x
x_vals = np.arange(1, n_genomes + 1)

# Plotting
plt.figure(figsize=(9, 6))

# Boxplots for pangenome
plt.boxplot([pan_all[:, i] for i in range(n_genomes)], positions=x_vals - 0.15, widths=0.25,
            patch_artist=True, boxprops=dict(facecolor='lightblue'), medianprops=dict(color='blue'), labels=['']*n_genomes)

# Boxplots for core genome
plt.boxplot([core_all[:, i] for i in range(n_genomes)], positions=x_vals + 0.15, widths=0.25,
            patch_artist=True, boxprops=dict(facecolor='lightgreen'), medianprops=dict(color='green'), labels=['']*n_genomes)

# Plot mean lines
plt.plot(x_vals, np.mean(pan_all, axis=0), label='Mean Pan-gene', color='blue', marker='o')
plt.plot(x_vals, np.mean(core_all, axis=0), label='Mean Core-gene', color='green', marker='o')

# Labels and aesthetics
plt.xlabel('Number of Haplotypes', fontsize=12, fontname="Arial")
plt.ylabel('Number of Orthogroups', fontsize=12, fontname="Arial")
plt.title('Pangenome and Core Genome Boxplots over 100 Permutations')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.xticks(x_vals, x_vals)
plt.tight_layout()
plt.savefig("pangenome_coregenome_boxplot.pdf", format="pdf")
plt.show()

