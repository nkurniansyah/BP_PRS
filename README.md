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

1.  Trait: SBP and DBP  
2.  Ancestry: GWAS Ancestry (which ancestry the GWAS summary statistics
    matched)  
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

We also provide the weight that we built in the Mass General Brigham
biobank to construct the PRS. The table below provides, for each
ancestry, trait summary statistic, and weight used for race/ethnicity,
the following information:  
1. Trait: SBP and DBP  
2. Ancestry: GWAS Ancestry (which ancestry the GWAS summary statistics
matched)  
3. All: Weight for all the sample  
4. Asian: Weight for the Asian population  
5. Balck: Weight for the Black population  
6. Hispanic/Latino: Weight for The Hispanic/Latino population  
7. White: Weight for the White population  

    ##   Trait Ancestry  All Asian Black Hispanic/Latino White
    ## 1   SBP      AFR 2.23 -0.01  3.92            3.54  2.15
    ## 2            EAS 0.84  1.60  0.55            0.76  0.83
    ## 3            EUR 4.03  3.58  3.11            4.37  4.06
    ## 4   DBP      AFR 0.97  0.61  2.30            1.31  0.85
    ## 5            EAS 0.49  1.13  0.46            0.76  0.47
    ## 6            EUR 2.13  2.67  1.78            2.48  2.11

## PLINK command for PRS construction

This command is to construct PRS using summary statistics that we
provide. The summary statistics are already based on the PRS-CSx
results. Note that genetic data files need to be specified in the –bfile
argument.

    for ancestry in {AFR,EAS,EUR}; do plink --bfile genetic/file/here \
    --score ../Summary_stats/SBP_$ancestry\_hg38.txt \
    --out  ../PRS/PRS_SBP_$ancestry ; done

## Constructing PRSsum based on trait and ancestry specific PRS

After constructing trait and ancestry-specific PRS, the BP-PRS is
obtained via the weighted PRSsum approach: First, we scaled PRS using
TOPMed mean and SD value of each trait and ancestry-specific PRS;
Second, we applied the weight from the MGB biobank that we provide for
each trait and ancestry-specific PRS; Third, we sum the PRS for each
ancestry for a specific trait. Finally, we applied the final scale for
BP-PRS using the TOPMed mean and SD we provided below.  

    ##    Trait            Race      Mean   SD
    ## 1    SBP             All -2.85e-14 3.22
    ## 2    SBP           Black -2.18e-14 3.05
    ## 3    SBP           Asian -2.65e-14 4.00
    ## 4    SBP Hispanic/Latino -3.07e-14 3.45
    ## 5    SBP           White -2.87e-14 3.25
    ## 6    DBP             All  5.92e-15 1.95
    ## 7    DBP           Black  7.61e-15 2.26
    ## 8    DBP           Asian  6.98e-15 2.66
    ## 9    DBP Hispanic/Latino  7.44e-15 2.34
    ## 10   DBP           White  5.65e-15 1.93
