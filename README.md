## Introduction

This repository provides information regarding the construction of a
polygenic risk score (PRS) for Systolic and Diastolic Blood Pressure
that we developed (BP-PRS) in the manuscript Evaluating the use of blood
pressure polygenic risk scores based on the largest available GWAS
across race/ethnic background groups (link to be added)

First, it provides instructions for constructing the BP-PRS based on
weighted summary statistics from PRS-CSx. These files can be downloaded
from the PGS catalog \[link to be added\], and code for using them to
construct the PRS. Second, this repository also provides code that we
used for the analyses in the manuscript (see folder “Code”)

## Required packages

We used [PLINK v1.9](https://www.cog-genomics.org/plink/ "PLINK v1.9")
to generate PRS. We provide example code that also uses to construct the
PRS.

    install.packages("dplyr")

Other software and packages that we used, but may not be necessary for
others to construct the PRS, are as follows:  
1. We performed the analysis using R version 4.0.2.  
2. We used the following packages from CRAN: dplyr, tidyverse,
data.table, purrr, pROC.  
3. We used the following packages from BioConductor: GENESIS,
GWASTools.  

## PRS construction

Our (SBP and DBP) BP-PRS is a weighted sum of multiple specifics
ancestry GWAS. we used and followed the
[PRS-CSx](https://github.com/getian107/PRScsx "PRS-CSx") to train
ancestry-specific effect sizes and create ancestry-specific PRS and
weighted sums of them. PRS-CSx takes summary statistics from multiple
GWAS. Each GWAS is assigned an ancestry and a reference panel that
matches this ancestry. Here, we paired
[UKBB+ICBP](https://www.nature.com/articles/s41588-018-0205-x,%22UKBB+ICBP%22)
with European ancestry,
[BBJ](https://www.nature.com/articles/s41588-018-0047-6 "BBJ") with East
Asian, and
[COGENT](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006728 "COGENT")
with the African reference panel. We used UKBB as LD reference
panels.See manuscript for more detail. Additionally, the code, TOPMed
mean, sd and weight from MGB Biobank to construct the PRS score are
available in this repository.  

The table below provides, for each ancestry and trait summary statistic
used, the following information:  

1.  Ancestry: GWAS Ancestry (which cohort/study the GWAS summary
    statistics matched)  
2.  Trait: (SBP, DBP)  
3.  TOPMed\_mean: the mean of the PRS after it was constructed in the
    multi-ethnic TOPMed population. That is, each of the TOPMed
    participants had a PRS value. This is the mean of these values  
4.  TOPMed\_sd: the standard deviation (SD) of the PRS after it was
    constructed in the multi-ethnic TOPMed population. That is, each of
    the TOPMed participants had a PRS value. This is the SD of these
    values  

<!-- -->

    ##   Trait Ancestry TOPMed_mean TOPMed_sd
    ## 1   SBP      EAS    4.20e-07  2.43e-07
    ## 2   SBP      AFR    9.29e-07  2.38e-06
    ## 3   SBP      EUR    1.49e-06  2.55e-07
    ## 4   DBP      EAS    2.04e-07  2.46e-07
    ## 5   DBP      AFR   -1.49e-06  1.94e-06
    ## 6   DBP      EUR    7.17e-07  2.19e-07
