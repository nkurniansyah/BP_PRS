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

1.  Trait: SBP and DBP  
2.  Ancestry: GWAS Ancestry (which ancestry the GWAS summary statistics
    matched)  
3.  All: Weight for all the sample  
4.  Asian: Weight for the Asian population  
5.  Balck: Weight for the Black population  
6.  Hispanic/Latino: Weight for The Hispanic/Latino population  
7.  White: Weight for the White population  

<!-- -->

    ##   Trait Ancestry  All Black Hispanic/Latino Asian White
    ## 1   SBP      AFR 1.81  2.00            1.78  0.48  1.82
    ## 2   SBP      EAS 0.68  0.35            0.50  1.17  0.69
    ## 3   SBP      EUR 3.28  2.31            3.37  2.58  3.35
    ## 4   DBP      AFR 0.69  1.42            0.76  0.63  0.61
    ## 5   DBP      EAS 0.40  0.38            0.53  0.83  0.39
    ## 6   DBP      EUR 1.71  1.44            1.82  1.93  1.71

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
    ## 1    SBP             All -2.32e-14 2.62
    ## 2    SBP           Black -1.62e-14 1.84
    ## 3    SBP           Asian -1.91e-14 2.64
    ## 4    SBP Hispanic/Latino -2.35e-14 2.63
    ## 5    SBP           White -2.36e-14 2.67
    ## 6    DBP             All  4.60e-15 1.57
    ## 7    DBP           Black  5.39e-15 1.56
    ## 8    DBP           Asian  5.41e-15 1.90
    ## 9    DBP Hispanic/Latino  5.07e-15 1.70
    ## 10   DBP           White  4.45e-15 1.57

See code below to construct weighted PRSsum.

    source("./Code/construct_wPRSsum.R")

    races<-c("All","Black","Asian","Hispanic/Latino","White")
    for(race in races){
      traits<-c("SBP","DBP")
      for(trait in traits){
        PRS_AFR<-paste0("../PRS/",trait,"_AFR.all_score")
        PRS_EAS<-paste0("../PRS/",trait,"_EAS.all_score")
        PRS_EUR<-paste0("../PRS/",trait,"_EUR.all_score")

        wprssum_score<-construct_wPRSsum(ethnic_background = race ,
                       TOPMed_scaling = TOPMed_scaling,
                       TOPMed_PRSsum_scaling = TOPMed_PRSsum_scaling,
                       mgb_weight=mgb_weight,
                       PRS_AFR_file=PRS_AFR,
                       PRS_EAS_file=PRS_EAS,
                       PRS_EUR_file=PRS_EUR)
        write.csv(wprssum_score, file = paste0("../output/wPRSsum/here/",trait,"_wPRSsum_",race,".csv"), row.names = F)
      }
      
    }

## Example code for association analysis

We performed association analysis using mixed models implemented in the
GENESIS R package for all individuals and linear regression for
unrelated individual. Below is an example code. It uses function that we
provide in the folder “Code”.

    library(GENESIS)
    library(GWASTools)
    library(pROC)


    source("./Code/*")


    #phenotype for all individual, if you want to run by background, please subset phenotype and PRSsum based on race/ethnicity.

    pheno<- fread(phenotype_file, data.table=F)


    # merge PRSsum with phenotype

    pheno_df<-left_join(pheno,prssum, by="person_id")


    covarites_prs<- c("BMI","age","sex","site","race",paste0("PC_",1:11),"wPRSsum")

    outcome<-"SBP"

    ## Kinship matrix

    covMatlist<-getobj(covMatlist)


    assoc_df<- run_assoc_mixmodel(pheno=pheno,
                                  outcome=outcome,
                                  covars_prs=covarites_prs, 
                                  covmat=covMatlist,
                                  group.var=NULL)


    # Perform AUC

    #only use unrelated people

    unrels<- getobj(unrels_people)

    pheno_unrels<- pheno[pheno_df$sample.id %in% unrels,]


    pve<-boot_pve(phenotype=pheno_unrels,
                  covariates_string=covarites_prs,
                  outcome=outcome,exposure="wPRSsum",seed=NULL, n=1000)




    final_assoc<-c(assoc_df,pve)
