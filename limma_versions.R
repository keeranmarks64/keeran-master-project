### Script to apply specific limma version to feature and label matrices

apply_limma_version <- function(df.test,label,limma_version){
  #' Function to apply a specific version of limma procedure to data and label matrices.
  #' 
  #' 
  #' @param df.test: data matrix of count values for each species and sample
  #' @param label: vector of length equivalent to the number of samples, representing the samples that are part of the case (1) and control (-1) groups. 
  #' @param limma_version: a string specifiying the specific version of the limma procedure to be applied. 
  #' @return returns the fitted limma object and p-value results matrix for the data provided. 
  
  if(limma_version == "original"){
    
    fit = limma::lmFit(df.test,design=data.frame(intercept = -1, label))
    res = limma::eBayes(fit)
    
  } else {
    
    stop("limma version not recognised")
  }
  
  result=list(fit,res)
  names(result)=c("fit","res")
  
  return(result)
  
}