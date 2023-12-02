# Sabolic_et_al_2023
As set of custom Perl and R scripts used for data analysis in Sabolic et al. (2023): Plastic and genomic change of a newly established lizard population following a founder event.

## coverage_vcfilter.pl
Perl script to analyse mean coverage/depth (DP) per locus.

## fix_fastq_headers.pl
Pearl script to fix headers in fastq files.

## genomic_PCA.R
R script for running PCA analysis on genomic datasets.

## HWeq_filter.R
R script to find and output a list of loci not in HWe in a given dataset.

## MCMCglmm_GeneticGroups.R
R script to calculate heritability using Bayesian MCMCglmm animal model with genetic groups.

## pop-pairwise_Fst.R
R script to estimate population-pairwise Fst values.

## run_diveRsity.R
R script to estimate genomic diversity indices.

## run_PCAdapt.R
R script for running PCAdapt genome scan to detect loci under selection.

## RunAndPlot_LFMM.R
R script for for running LFMM analysis and extracting and plotting the results.

## RunAndPlot_RDA.R
R script for for running RDA analysis on genomic data and extracting and plotting the results.

## standardize_fastq_length.pl
Perl script to standardize length of sequences in fastq files.

## vcf_random_resampler.pl
Perl script to randomly subsample VCF file to N samples per population.

## vcfixer.pl
Perl script to check VCF file, remove loci and samples with too many missing, and input remaining missing genotypes.
