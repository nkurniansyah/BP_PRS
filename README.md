## Introduction

This repository provides information regarding the construction of a
polygenic risk score (PRS) for Systolic and Diastolic Blood Pressure
that we developed (BP-PRS) in manuscript Evaluating the use of blood
pressure polygenic risk scores based on largest available GWAS across
race/ethnic background groups (link to be added)

First, we used and followed the
[PRS-CSx](https://github.com/getian107/PRScsx "PRS-CSx") to train
ancestry-specific effect sizes and create ancestry-specific PRS, and
weighted sums them. PRS-CSx takes summary statistics from multiple GWAS.
Each GWAS is assigned an ancestry and a reference panel that matches
this ancestry. Here, we paired UKBB+ICBP with the European ancestry
reference panel, and BBJ with the East Asian panel, and MVP with the
African reference panel. Additionally, the code, and weight from MGB
Biobank to construct the PRS score are availbale in this repository.

Second, this repository also provides code the we used for the analyses
in the manuscript (see folder “Code”).

## Required packages

We used [PLINK v1.9](https://www.cog-genomics.org/plink/ "PLINK v1.9")
to generate PRS. We provide example code that also uses to construct the
PRS.

Other software and packages that we used, but may not be necessary for
others to construct the PRS, are as follows:  
1. We performed the analysis using R version 4.0.2.  
2. We used the following packages from CRAN: dplyr, tidyverse,
data.table, purrr, pROC.  
3. We used the following packages from BioConductor: GENESIS,
GWASTools.  

    install.packages("dplyr")

## PRS construction

Our BP-PRS is a weighted sum of multiple specifics ancestry GWAS.
Summary statistics to create the each of the ancstry -PRS are provided
in here.

The table below provides, for each ancestry-specific GWAS used, the
following information:  

1.  GWAS\_Anc: GWAS Ancestry (which cohort/study the GWAS summary
    statistics matched)  
2.  Trait (SBP, DBP)  
3.  TOPMed\_mean: the mean of the PRS after it was constructed in the
    multi-ethnic TOPMed population. That is, each of the TOPMed
    participants had a PRS value. This is the mean of these values.  
4.  TOPMed\_sd: the standard deviation (SD) of the PRS after it was
    constructed in the multi-ethnic TOPMed population. That is, each of
    the TOPMed participants had a PRS value. This is the SD of these
    values.

<!-- -->

    ##   GWAS_Anc Trait TOPMed_mean TOPMed_sd
    ## 1      AFR   SBP   -6.73e-07  4.54e-07
    ## 2      EUR   SBP    1.31e-06  3.22e-07
    ## 3      EAS   SBP    4.99e-07  3.05e-07
    ## 4      AFR   DBP   -4.02e-07  3.05e-07
    ## 5      EUR   DBP    6.64e-07  2.60e-07
    ## 6      EAS   DBP    3.02e-07  2.45e-07

We also provide MGB Biobank-trained PRS summation weights for the best
performing PRS.

The table below provides, for each background-specific weight used to
construct PRS, the following information:  

1.  Trait (SBP, DBP)  
2.  Race\_Ethnic\_Background: Race/ ethnic background
3.  MBG\_weight: weight in MGB to construct PRS.

<!-- -->

    ##   GWAS_Anc Trait TOPMed_mean TOPMed_sd
    ## 1      AFR   SBP   -6.73e-07  4.54e-07
    ## 2      EUR   SBP    1.31e-06  3.22e-07
    ## 3      EAS   SBP    4.99e-07  3.05e-07
    ## 4      AFR   DBP   -4.02e-07  3.05e-07
    ## 5      EUR   DBP    6.64e-07  2.60e-07
    ## 6      EAS   DBP    3.02e-07  2.45e-07

## PRSice command for PRS construction

This command is to construct PRS using the summary statistics that we
provide. No clumping is needed and no selection of SNPs. The summary
statistics are already based on the specific set of SNPs selected after
clumping and setting a p-value threshold. Note that genetic data files
need to be specified in the –target argument.



    Rscript ./PRSice.R \
     --dir ./PRS_Output \
     --prsice ./PRSice_linux/PRSice_linux \
     --base ./Summary_Statistics_for_PRS_construction/. \
     --target ./Genotype \
     --thread 2 \
     --chr Chromosome 
     --bp Position 
     --A1 Allele1 
     --A2 Allele2 
     --pvalue PValue \
     --bar-levels Threshold \
     --stat BETA 
     --all-score T \
     --out ./out_prs \
     --no-clump T
     --print-snp T \
     --ignore-fid T 
     --no-regress T 
     --fastscore T 
     --model add 
     --no-full T 
     --chr-id c:l:a:b

## Constructing PRSsum based on trait-specific PRS

After constructing trait-specific PRS, the HTN-PRS is obtained via the
PRSsum approach: as an unweighted of the scaled trait-specific PRS. For
scaling, we use the TOPMed mean and SD values of each trait-specific
PRS, and we also provide here the TOPMed mean and SD of the HTN-PRS for
final scaling. Using the same scaling throughout guarantees that effect
size estimates are similarly interpreted across all datasets and
individuals who use this PRS.

    ##   TOPMed_mean TOPMed_sd
    ## 1   -4.96e-16      2.69

See code below to construct PRSsum.

    library(data.table)
    library(dplyr)
    library(purrr)


    prs_traits <- c("SBP", "DBP","HTN")
    out<-list()
    for(trait in prs_traits){
      
      
      prs_output <-paste0("./", prs,".txt")
      prs_df <-fread(prs_output, data.table=F)
      prs_df <- prs_df %>% dplyr::select(-IID)
      colnames(prs_df)<- c("sample.id", trait)
      
      #standardize trait-prs using the mean and sd from TOPMed 
      cur_TOPMed_mean <- PRS_info$TOPMed_mean[which(PRS_info$Trait == trait)]
      cur_TOPMed_sd <- PRS_info$TOPMed_sd[which(PRS_info$Trait == trait)]
      
      prs_df[, trait]<- (prs_df[, trait] - cur_TOPMed_mean)/cur_TOPMed_sd
      out[[trait]]<- prs_df
      
      
    }

    combine_prs <- purrr::reduce(out, full_join , by="sample.id")



    prssum<- data.frame(sample.id=prssum$sample.id, 
                        PRSsum=apply(prssum[,], 1, sum))

    prssum[,"PRSsum"]<- (prssum[,"PRSsum"] - TOPMed_HTN_PRS_mean_sd$TOPMed_mean))/TOPMed_HTN_PRS_mean_sd$TOPMed_sd

## Example code for association analsis

We performed association analysis using mixed models implemented in the
GENESIS R package. Below is an example code. It uses function that we
provide in the folder “Code”.

    library(GENESIS)
    library(GWASTools)
    library(pROC)


    source("./Code/*")


    #phenotype

    pheno<- fread(phenotype_file, data.table=F)


    # merge PRSsum with phenotype

    pheno_df<-left_join(pheno,prssum, by="sample.id" )


    covarites_prs<- c("BMI","age","sex","site","race",paste0("PC_",1:11),"PRSsum")

    outcome<-"HTN"

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


    auc<- generate_auc(pheno=pheno_unrels,
                       outcome=outcome,
                       covars_prs=covarites_prs, seed=NULL,
                       n= 2000)



    final_assoc<-c(assoc_df,auc)
