#' Title : Generate Bootstrap PVE
#'
#' @param pheno : data frame of the phenotype (Unrelated individuals)
#' @param outcome : as.numeric , outcome to test 
#' @param covars_prs : covariates to adjust 
#' @param n : boostrap n
#' @param seed : random seed number
#'
#' @return CI pve
#' @export
#'
#' @examples
#' 

boot_pve <- function(phenotype,covariates_string,outcome,exposure,seed=NULL, n=1000){
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(!any(grepl("prs_", exposure))){
    exp<-paste0(exposure)
    #xposure<-paste0("`",exposure,"`")
    exposure<-exposure
  }else{
    exposure<-exposure
  }
  
  message(paste("using ", exposure))
  covars_no_prs<- covariates_string
  covars_prs<- c(covars_no_prs,exposure )
  
  
  head(phenotype)
  #print("yes")
  colnames(phenotype)
  phenotype_clean<- phenotype[,c("sample.id",outcome,covars_no_prs,exposure)]
  #phenotype_clean<- phenotype %>% dplyr::select(sample.id,outcome,covars_no_prs,exposure)
  #print("no")
  head(phenotype_clean)
  dim(phenotype_clean)
  phenotype_clean[,exposure]<- sapply(phenotype_clean[,exposure], as.numeric)
  
  #phenotype_clean[,exposure][is.nan(phenotype_clean[,exposure])]<-0
  
  ## Checking binary
  phenotype_clean<-na.omit(phenotype_clean) 
  head(phenotype_clean)
  dim(phenotype_clean)
  # print(covars_no_prs)
  # print(covars_prs)
  
  
  out <- list()
  
  print("here")
  for(i in 1:n){
    
    df <- sample_n(phenotype_clean, nrow(phenotype_clean), replace = TRUE)
    head(df)
    
    #print(paste0(outcome," is Gaussian..."))
    fit_no_prs<- (lm(as.formula(paste(outcome, "~", paste(covars_no_prs, collapse = "+" ))),data=df))
    
    fit_prs<- (lm(as.formula(paste(outcome, "~", paste(covars_prs, collapse = "+" ))),data=df))
    resid1<-mean(fit_no_prs$residuals^2)
    
    #fit_prs<- (lm(as.formula(paste(outcome, "~", paste(covars_prs, collapse = "+" ))),data=phenotype_clean))
    resid2<-mean(fit_prs$residuals^2)
    
    #print(covars_prs)
    variance_explained<- (resid1-resid2)/resid1 * 100
    out[[i]]<-variance_explained
  }
  
  res<-do.call(c,out)
  

  variance_explained<-quantile(res, c(0.025,0.5, 0.975))

  message("boostrap completed")
  return(variance_explained)
}