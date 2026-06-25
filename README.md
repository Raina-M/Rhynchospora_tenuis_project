# Rhynchospora_tenuis_project
Cutstom scripts for the manuscript of _R. tenuis_ project.

> Zhang, Meng, et al. "Sex without crossovers mimics clonal reproduction in the holocentric plant Rhynchospora tenuis." bioRxiv (2026): 2026-01.
> doi: https://doi.org/10.64898/2026.01.17.700054

## Syteny analyses
`syri_plot.sh`: DNA collinearity by [SyRI](https://github.com/schneebergerlab/syri) and [plotsr](https://github.com/schneebergerlab/plotsr). Reference sctipts for __Figure 1b__ and __Extended Figure 2__.

`run_genespace.R`: Genome syteny based on orthologous genes. Reference script for running [GeneSpace](https://github.com/jtlovell/GENESPACE) in __Figure 1a__.

## Divergence time estimation of two haplotypes
`divergence_time_estimation/Fig_1f_estimate_divergence_time.R`: Divergence time calculation and plot for __Figure 1f__.

- The input Ks values between two haplotypes was computed by [CoGe](https://ghibli.bti.cornell.edu/coge/) based on the CDS sequences annotated by
Helixer. All Ks of the CDSs overlapped with TE annotations were removed. Refer to the script `divergence_time_estimation/remove_TE_CDS.sh`
- The SNP input was derived from SyRI from above synteny analysis, then the SNP density (SNP per Mb) were calculated as the input.
- The LTR divergence used the output from SubPhaser. Please check `divergence_time_estimation/extract_LTR_insertion_time_from_SubPhaser.sh`

## Single cell analyses of pollen nuclei
The genotyping of single pollen nuclei used the same pipeline as described in [Castellani, M., Zhang, M. et al](https://doi.org/10.1038/s41477-024-01625-y). The plots in __Figure 3a__ was also part of its pipeline (`hapCO_identification.R` in the [pipeline](https://github.com/Raina-M/detectCO_by_scRNAseq) GitHub page).

`segregation_bias_ChiSquareTest.R`: Non-mendelian segregation validation in __Figure 3d__.

## Single cell ATAC-seq analyses of embryo and endosperm nuclei
The pipeline is similar to single cell pipeline for pollen nuclei in general with minor modifications. The entire pipeline is available in the folder `embryo_scATAC_analyses/`. This pipeline took the embryo _R. tenuis_ REC accession as an example, but the endosperm analyses has no difference to it. 

## F1 genotyping

## ChIP-seq pipeline and relevant analyses
`ChIP_and_methylation_metaplots/` has the pipeline for DNA methylation data analyses using [Bismark](https://github.com/FelixKrueger/Bismark) and other following plotting scripts. The other ChIP-seq data and analyses were in our another publication [Hofstatter et al.](https://www.sciencedirect.com/science/article/pii/S0092867422007978).

`2_metaplot.sh`: Metaplots in __Extended Data Fig. 1c,d__.

## Phylogenetic relation of _R. tenuis_ accessions 

## Pseudogene analyses

## Gene conversion detection based on F1 WGS data

## Other
__Extended Data Fig. 4a__
__Extended Data Fig. 5a,b__
__Supplementary Fig. 14b, d__


## Final comment
 please open an issue if you have any doubts about the scripts and data processing questions regarding our publication.
