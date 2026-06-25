# Rhynchospora_tenuis_project
Cutstom scripts for the manuscript of _R. tenuis_ project.

> Zhang, Meng, et al. "Sex without crossovers mimics clonal reproduction in the holocentric plant Rhynchospora tenuis." bioRxiv (2026): 2026-01.
> doi: https://doi.org/10.64898/2026.01.17.700054

### Syteny analyses
`syri_plot.sh`: DNA collinearity by [SyRI](https://github.com/schneebergerlab/syri) and [plotsr](https://github.com/schneebergerlab/plotsr). Reference sctipts for __Figure 1b__ and __Extended Figure 2__.

`run_genespace.R`: Genome syteny based on orthologous genes. Reference script for __Figure 1a__.

### Divergence time estimation of two haplotypes
`divergence_time_estimation/Fig_1f_estimate_divergence_time.R`: Divergence time calculation and plot for __Figure 1f__.

- The input Ks values between two haplotypes was computed by [CoGe](https://ghibli.bti.cornell.edu/coge/) based on the CDS sequences annotated by
Helixer. All Ks of the CDSs overlapped with TE annotations were removed. Refer to the script `divergence_time_estimation/remove_TE_CDS.sh`
- The SNP input was derived from SyRI from above synteny analysis, then the SNP density (SNP per Mb) were calculated as the input.
- The LTR divergence used the output from SubPhaser. Please check `divergence_time_estimation/extract_LTR_insertion_time_from_SubPhaser.sh`

