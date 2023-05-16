

#' Title: Construct weighted PRSsum
#'
#' @param ethnic_background: race/ethnicity (All, Black, White,Hispnic/Latino, Asian)
#' @param TOPMed_scaling : TOPMed mean and SD in data frame for each Trait and PRS ancestry
#' @param trait : BP trait selection (SBP and DBP)
#' @param mgb_weight : MGB weight in data frame for each trait, ancestry and race/ethnicity
#' @param TOPMed_PRSsum_scaling : Final scaling for PRSsum 
#' @param PRS_AFR_file : PRS AFR file 
#' @param PRS_EAS_file : PRS EAS file 
#' @param PRS_EUR_file : PRS EUR file 
#'
#' @return PRSsum for selection trait and race/ethnicity
#' @export
#'
#' @examples
#' 
construct_wPRSsum<- function(ethnic_background, TOPMed_scaling,
                             trait, mgb_weight, 
                             TOPMed_PRSsum_scaling,
                             PRS_AFR_file,
                             PRS_EAS_file,
                             PRS_EUR_file){
  
  
  TOPMed_scaling_selected<-TOPMed_scaling[which(TOPMed_scaling$Trait==trait),]
  mgb_weight_selected<-mgb_weight[which(mgb_weight$Trait==trait),c("Trait","Ancestry",ethnic_background)]
  colnames(mgb_weight_selected)<-c("Trait","Ancestry","MGB_Weight")
  
  
  #ancestries<-unique(as.character(TOPMed_scaling_selected$Ancestry))
  
  prs_file<-list(AFR=PRS_AFR_file, EUR=PRS_EUR_file, EAS=PRS_EAS_file)
  
  for(ancestry in names(prs_file)){
    prs_ancestry<-prs_file[[ancestry]]
    mgb_weigt_ancestry<-mgb_weight_selected[which(mgb_weight_selected$Ancestry==ancestry),]
    
    prs_df<- fread(prs_ancestry, data.table = F)
    
    prs_df<-prs_df %>% dplyr::select(IID,SCORE)
    colnames(prs_df)<-c("person_id",paste0("prs_", ancesestry))
    
    prs_df[,2]<-as.numeric(prs_df[,2])
    
    topmed_scaling_ancestry<-TOPMed_scaling_selected[which(TOPMed_scaling_selected$Ancestry==ancesestry),]
    mean_topmed<- as.numeric(topmed_scaling_ancestry$TOPMed_mean)
    sd_topmed<- as.numeric(topmed_scaling_ancestry$TOPMed_sd)
    prs_df[,paste0("prs_",ancesestry)]<-(prs_df[,paste0("prs_",ancesestry)]-mean_topmed)/sd_topmed
    
    prs_df[,paste0("prs_",ancesestry)]<-prs_df[,paste0("prs_",ancesestry)]*mgb_weigt_ancestry$MGB_Weight
    
    out[[ancesestry]]<- prs_df
    
  }
  
 
  prs_all<-purrr:::reduce(out, left_join, by="person_id")
  
  prs_sum<-data.frame(person_id=prs_all$person_id, prs=apply(prs_all[,2:ncol(prs_all)],1, sum))
  TOPMed_PRSsum_scaling_selected<-TOPMed_PRSsum_scaling[which(TOPMed_PRSsum_scaling$Trait==trait & TOPMed_PRSsum_scaling$Self_reported_ethnic_background == ethnic_background),]

  prs_sum[,"prs"]<-(prs_sum[,"prs"]-TOPMed_PRSsum_scaling_selected$Mean)/TOPMed_PRSsum_scaling_selected$SD
  return(prs_sum)
}
