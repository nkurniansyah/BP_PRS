## Introduction

This repository provides information regarding the construction of a
polygenic risk score (PRS) for Systolic and Diastolic Blood Pressure
that we developed (BP-PRS) in the manuscript Evaluating the use of blood
pressure polygenic risk scores based on the largest available GWAS
across race/ethnic background groups (link to be added) First, we used
and followed the[PRS-CSx](https://github.com/getian107/PRScsx "PRS-CSx")
to train ancestry-specific effect sizes and create ancestry-specific PRS
and weighted sums of them. PRS-CSx takes summary statistics from
multiple GWAS. Each GWAS is assigned an ancestry and a reference panel
that matches this ancestry. Here, we paired UKBB+ICBP with European
ancestry, BBJ with East Asian, and COGENT with the African reference
panel. Additionally, the code and weight from MGB Biobank to construct
the PRS score are available in this repository. Second, this repository
also provides code that we used for the analyses in the manuscript (see
folder “Code”)

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

\#\#PRS construction Our BP-PRS is a weighted sum of multiple specifics
ancestry GWAS. Summary statistics to create the each of the ancstry -PRS
are provided in PGS catalog \[link to be added\]

We provide the command to construct PRS based summary statistics that we
provide. Genetic data files need to be specified in the –bfile argument.
Below is an example command to construct PRS for SBP.

    #
    for ancestry in {AFR,EAS,EUR}; do plink --bfile genetic/file/here \
        --score ../Summary_stats/SBP_$ancestry\_hg38.txt \
        --out  ../PRS/PRS_SBP_$ancestry ; done

The table below provides, for each ancestry and trait summary statistic
used, the following information:

Ancestry: GWAS Ancestry (which cohort/study the GWAS summary statistics
matched)  
Trait (SBP, DBP)  
Threshold: p-value threshold for selecting SNPs into the PRS
TOPMed\_mean: the mean of the PRS after it was constructed in the
multi-ethnic TOPMed population. That is, each of the TOPMed participants
had a PRS value. This is the mean of these values. TOPMed\_sd: the
standard deviation (SD) of the PRS after it was constructed in the
multi-ethnic TOPMed population. That is, each of the TOPMed participants
had a PRS value. This is the SD of these values.

    ##   Trait Ancestry      Mean       SD
    ## 1   SBP      EAS  4.20e-07 2.43e-07
    ## 2            AFR  9.29e-07 2.38e-06
    ## 3            EUR  1.49e-06 2.55e-07
    ## 4   DBP      EAS  2.04e-07 2.46e-07
    ## 5            AFR -1.49e-06 1.94e-06
    ## 6            EUR  7.17e-07 2.19e-07

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

## Constructing weighted PRS sum based on ancestry-specific PRS

After constructing trait-specific PRS, the BP-PRS is obtained via the
weigted PRS sum approach: as an weighted of the scaled ancestry-specific
PRS. For scaling, we use the TOPMed mean and SD values of each
ancestry-specific PRS, and we also provide here the TOPMed mean and SD
of the BP-PRS for final scaling for each ancestry. Using the same
scaling throughout guarantees that effect size estimates are similarly
interpreted across all datasets and individuals who use this PRS.

    ##    Trait     Race_Ethnic TOPMed_mean TOPMed_sd
    ## 1    SBP   Multi-Ethninc   -4.88e-16      3.01
    ## 2    SBP           White   -5.02e-16      3.12
    ## 3    SBP           Black   -2.68e-16      1.83
    ## 4    SBP           Asian   -3.78e-16      1.82
    ## 5    SBP Hispanic/Latino   -4.65e-16      2.78
    ## 6    DBP   Multi-Ethninc   -3.67e-16      1.72
    ## 7    DBP           White   -3.72e-16      1.74
    ## 8    DBP           Black   -2.51e-16      1.27
    ## 9    DBP           Asian   -3.61e-16      1.61
    ## 10   DBP Hispanic/Latino   -3.96e-16      1.77

See code below to construct PRSsum.

    library(data.table)
    library(dplyr)
    library(purrr)

    create_prs<- function(ethnic_background, PRS_info,trait, mgb_weight_files){
        

        PRS_info_trait<-PRS_info[which(PRS_info$Trait==trait),]
      
        ancestries<-unique(as.character(topmed_mean_sd_file$GWAS_Anc))

        out<-list()
        for(ancesestry in ancestries){
        
            prs_file<- paste0("../PRS/PRS_",trait,"_"ancestry,".profile")
            prs_df<- fread(prs_file, data.table = F)

            prs_df<-prs_df %>% dplyr::select(IID,SCORE)
            colnames(prs_df)<-c("person_id",paste0("prs_", ancesestry))

            prs_df[,2]<-as.numeric(prs_df[,2])
            
            info<-PRS_info_trait[which(PRS_info_trait$GWAS_Anc==ancesestry),]
            mean_topmed<- as.numeric(info$TOPMed_mean)
            sd_topmed<- as.numeric(info$TOPMed_sd)
            prs_df[,paste0("prs_",ancesestry)]<-(prs_df[,paste0("prs_",ancesestry)]-mean_topmed)/sd_topmed
            out[[ancesestry]]<- prs_df
        }
        prs_all<-reduce(out, left_join, by="person_id")
        #mgb_weight_file<-paste0("/Users/nuzululkurniansyah/Documents/BP_PRS/MGB_Weight/PRS-CSx1_Weight.csv")
        mgb_weight<-fread(mgb_weight_file, data.table = F)
        mgb_weight<-mgb_weight[which(mgb_weight$Trait==trait),]
        mgb_weight<-mgb_weight[,c("Trait","Ancestry",ethnic)]
        out_weight<-list()
        for(i in 1:nrow(mgb_weight)){
            mgb_weight_df<- mgb_weight[i,]
            prs_name<-paste0("prs_",mgb_weight_df$Ancestry)
            prs_val<-prs_all[,c("person_id",prs_name])]
            prs_val[,prs_name]<-prs_val_df[,prs_name]*as.numeric(mgb_weight_df[,ethnic])
            out_eff[[i]]<-prs_val
        }
        prs_weight<-reduce(out_eff, left_join, by="person_id")
        
        prs_sum<-data.frame(person_id=prs_weight$person_id, prs=apply(prs_weight[,2:ncol(prs_weight)],1, sum))

        prs_sum_info<-topmed_mean[,c("type","prs")]
        prs_sum_info
        mean_topmed_prsum<- as.numeric(prs_sum_info[which(prs_sum_info$type=="mean_topmed"),2])
        sd_topmed_prsum<- as.numeric(prs_sum_info[which(prs_sum_info$type=="sd_topmed"),2])

        prs_sum_scale<-prs_sum
        prs_sum_scale[,"prs"]<-(prs_sum_scale[,"prs"]-mean_topmed_prsum)/sd_topmed_prsum
        head(prs_sum_scale)
        return(prs_sum_scale)
    }

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
