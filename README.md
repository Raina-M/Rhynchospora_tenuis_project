# Rhynchospora_tenuis_project
Cutstom scripts for the manuscript of _R. tenuis_ project.

> Zhang, Meng, et al. "Sex without crossovers mimics clonal reproduction in the holocentric plant Rhynchospora tenuis." bioRxiv (2026): 2026-01.
> doi: https://doi.org/10.64898/2026.01.17.700054

## Synteny analyses
`syri_plot.sh`: DNA collinearity by [SyRI](https://github.com/schneebergerlab/syri) and [plotsr](https://github.com/schneebergerlab/plotsr). Reference sctipts for __Figure 1b__ and __Extended Figure 2__.

`run_genespace.R`: Genome syteny based on orthologous genes. Reference script for running [GeneSpace](https://github.com/jtlovell/GENESPACE) in __Figure 1a__.


## Divergence time estimation of two haplotypes
`divergence_time_estimation/Fig_1f_estimate_divergence_time.R`: Divergence time calculation and plot for __Figure 1f__.

- The input Ks values between two haplotypes was computed by [CoGe](https://ghibli.bti.cornell.edu/coge/) based on the CDS sequences annotated by
Helixer. All Ks of the CDSs overlapped with TE annotations were removed. Refer to the script `divergence_time_estimation/remove_TE_CDS.sh`
- The SNP input was derived from SyRI from above synteny analysis, then the SNP density (SNP per Mb) were calculated as the input.
- The LTR divergence used the output from SubPhaser. Please check `divergence_time_estimation/extract_LTR_insertion_time_from_SubPhaser.sh`


## Single cell analyses of pollen nuclei
The genotyping of single pollen nuclei used the same pipeline as described in [Castellani, M., Zhang, M. et al](https://doi.org/10.1038/s41477-024-01625-y). The plots in __Figure 3a__ was also part of this pipeline (`hapCO_identification.R` in the [pipeline](https://github.com/Raina-M/detectCO_by_scRNAseq) GitHub page).

`segregation_bias_ChiSquareTest.R`: Non-mendelian segregation validation in __Figure 3d__.


## Single cell ATAC-seq analyses of embryo and endosperm nuclei
The pipeline is similar to single cell pipeline for pollen nuclei in general with minor modifications. The complete pipeline is available in the folder `embryo_scATAC_analyses/`. This pipeline took the 10X scATAC-seq of embryo nuclei from _R. tenuis_ REC accession as an example. The endosperm analyses has no difference to it. 

`plot_AF_along_chrs.R`: __Figure 4j__.

`plot_AF_along_chrs.R`: __Figure 4k__.


## F1 genotyping
The genotyping of selfed F1 plants used the same pipeline as described in [Castellani, M., Zhang, M. et al](https://doi.org/10.1038/s41477-024-01625-y). We uploaded an example run here with script name as `tiger_pipeline.sh`. Please note that the key input of this pipeline is ´VCF´ file that contains the SNPs of F1 WGS reads aligned to one of the haplotypes of _R. tenuis_ (REC). Input markers are genotyping markers defined by haplotype-specific SNPs, which was called from the alignments of deep-sequenced _R. tenuis_ (REC) short-reads to haplotype 1 assembly. 

`plottings/Fig_3b_plot_F1_AF.R`: Script of __Figure 3b__ for allele frequency along chromosomes and its density margin plots of a F1 individual.


## ChIP-seq pipeline and relevant analyses
`ChIP_and_methylation_metaplots/` includes the pipeline for DNA methylation data analyses using [Bismark](https://github.com/FelixKrueger/Bismark) and other following plotting scripts. The other ChIP-seq data and analyses were available in our previous publication [Hofstatter et al.](https://www.sciencedirect.com/science/article/pii/S0092867422007978).

`2_metaplot.sh`: Scripts for the metaplots in __Extended Data Fig. 1c,d__.

`4_correlation.R` Correlation analysis of epigentic signals in the metaplots (__Extended Data Fig. 1c,d__) - CENH3, H3K4me3, H3K9me2, CG, CHG, CHH between _R. pubera_, _R. breviuscula_, _R. tenuis_ (REC) in __Supplementary Table 2__.


## Phylogenetic relation of _R. tenuis_ accessions 

## Pseudogene analyses
`search_pseudogenes/` includes the complete pileline of pseudogene analyses with Helixer gene annotations and RNA alignments as inputs. RNA alignments were done by mapping RNA short-read sequencing data to diploid _R. breviuscula_ and _R. tenuis_ (REC) assembly using STAR. Please note that the results can differ a lot when including the CDSs of TEs. So it is highly recommended to remove the coding sequences of TEs when searching pseudogenes, which are usually largely accumulated in many plant genomes.

`03_report_pseudogenes.R`: Scripts for the plots in __Extended Data Fig. 5c__.


## Gene conversion detection based on F1 WGS data
`detect_gene_conversion.R`: The script shows how we detect gene conversion events. Input file is from our F1 genotyping pipeline with [TIGER](https://github.com/Imoteph/TIGER_Whole-Genome_Genotyping-by-Sequencing). The input file is `input_corrected.txt` from the first step of `tiger_pipeline.sh` in __#F1 genotyping__, which contains six fields: 
```
1 chrom  2 pos  3 hap1_allele  4 hap1_count  5 hap2_allele  6 hap2_count
```

Although the detection pipeline output some potential gene conversions, they are not reliable after careful manual check in IGV, because we never observed a single read that can support gene conversion events. The mean coverage of our F1 samples is 1x-3x. So our result does not exclude the possibility of gene conversion events in _R. tenuis_, considering the sequencing limitation of our data.


## Other
__Extended Data Fig. 5a,b__

__Supplementary Fig. 14b, d__


## Final comment
 please open an issue if you have any doubts about the scripts and data processing questions regarding our publication.
