# EGF-MDD-MR

## Summary
In this study, we performed two-sample Mendelian randomization (MR) analyses to examine the potential causal association between genetically predicted cortical expression of neurotrophic factor genes, including EGF, EGFR, FGF2, GDNF, and NGF, and major depressive disorder (MDD). Cis-eQTL instruments were obtained from the MetaBrain consortium, and GWAS summary statistics of MDD were obtained from the largest GWAS meta-analysis of depression to date (Als et al., 2023, Nat Med), with 23andMe data excluded due to access restrictions. <br>
Refer `00_Process.sh` for the major processes.

## Prepare the cis-eQTL data
1. Download the Cortex-EUR cis-eQTLs of the neurotrophic factor genes from the [MetaBrain consortium](https://www.metabrain.nl/cis-eqtls.html) (De Klein et al., 2023, Nat Genet).
2. Curate the cis-eQTLs of EGF, EGFR, FGF2, GDNF, and NGF into the instrumental variable table `data/mr_MetaBrain.csv` (columns: Symbol, SNP, A1, A2, CHR, BP, BETA, SE, P, EAF, Samplesize, F_statistics). BDNF was not analyzed as no cis-eQTL instruments were obtained. For more details, please refer to our Manuscript and sTable 2.

## Prepare the MDD GWAS data
1. GWAS summary statistics of MDD were obtained from Als et al., 2023, Nat Med, including 294,322 cases and 741,438 controls after excluding 23andMe data. The summary statistics can be requested from the [iPSYCH](https://ipsych.dk/en/research/downloads) website.
2. Place the summary statistics at `data/MDD.txt.gz` (tab-separated, columns: SNP, A1, A2, BETA, SE, P, N, CHR, BP).

## MR analysis
We conducted the MR analyses using the [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/index.html) package. Wald ratio and inverse variance weighting (IVW) were used as the primary methods, with MR-RAPS as an additional approach, and MR Egger, weighted median, and weighted mode for sensitivity analyses when sufficient instrumental variables were available. Genetic instruments were harmonized with default settings, and those with F-statistics < 10 were removed. Odds ratios with 95% confidence intervals were generated, and Benjamini-Hochberg correction was applied across tests.
1. Conduct the MR analysis using script `01_MR.R`
2. Prepare `data/plot_EGF_mr.res.csv` from the MR results with formatted labels and BH-adjusted q-values, and plot the forest plot using script `02_Forest.R`

## Literature mining
We explored the literature-based associations between EGF/EGFR and depression using [MELODI Presto](https://melodi-presto.mrcieu.ac.uk/). The results were shown in sTable 3 of our Manuscript.

## Dependencies
R packages: TwoSampleMR, mr.raps, ieugwasr, gwasvcf, RadialMR, tidyverse, data.table, plyr, openxlsx, R.utils, grid, forestploter, ragg. <br>
[PLINK](https://www.cog-genomics.org/plink/) with a 1000 Genomes EUR reference panel is only required if the LD proxy search block in `01_MR.R` is enabled.
