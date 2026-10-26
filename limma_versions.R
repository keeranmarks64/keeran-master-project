### Script to apply specific limma version to feature and label matrices

apply_limma_version <- function(df.test,label,limma_version){
  #' Function to apply a specific version of limma procedure to data and label matrices.
  #' 
  #' 
  #' @param df.test: data matrix of count values for each species and sample
  #' @param label: vector of length equivalent to the number of samples, representing the samples that are part of the case (1) and control (-1) groups. 
  #' @param limma_version: a string specifiying the specific version of the limma procedure to be applied. 
  #' @return returns the fitted limma object and p-value results matrix for the data provided. 
  
  if(limma_version == "limma_original"){
    
    #No Normalization
    fit = limma::lmFit(df.test,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit)
    
  } else if(limma_version == "limma_original_TMM"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying TMM Normalization
    dge = edgeR::normLibSizes(dge)
    
    #Converting counts to logCPM values
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=3)
    
    #Applying limma pipeline
    fit = limma::lmFit(logCPM,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit, trend=FALSE)
    
  } else if(limma_version == "limma_original_RLE"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying RLE Normalization
    dge = edgeR::normLibSizes(dge,method='RLE')
    
    #Converting counts to logCPM values
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=3)
  
    #Applying limma-trend pipeline
    fit = limma::lmFit(logCPM,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit, trend=FALSE)
    
  } else if(limma_version == "limma_original_TMMwsp"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying TMMwsp Normalization
    dge = edgeR::normLibSizes(dge,method='TMMwsp')
    
    #Converting counts to logCPM values
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=3)
    
    #Applying limma-trend pipeline
    fit = limma::lmFit(logCPM,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit, trend=FALSE)
    
  } else if(limma_version == "limma_trend"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying TMM Normalization
    dge = edgeR::normLibSizes(dge)
    
    #Converting counts to logCPM values
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=3)
    
    #Applying limma-trend pipeline
    fit = limma::lmFit(logCPM,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit, trend=TRUE)
    
  } else if(limma_version == "limma_trend_RLE"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying RlE Normalization
    dge = edgeR::normLibSizes(dge,method='RLE')
    
    #Converting counts to logCPM values
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=3)
    
    #Applying limma-trend pipeline
    fit = limma::lmFit(logCPM,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit, trend=TRUE)
    
  } else if(limma_version == "limma_trend_TMMwsp"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying TMMwsp Normalization
    dge = edgeR::normLibSizes(dge,method='TMMwsp')
    
    #Converting counts to logCPM values
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=3)
    
    #Applying limma-trend pipeline
    fit = limma::lmFit(logCPM,design=data.frame(intercept = 1, label=1-(label<0)))
    res = limma::eBayes(fit, trend=TRUE)
    
  } else if(limma_version == "limma_voom"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying TMM Normalization
    dge = edgeR::normLibSizes(dge)
    
    #Applying limma-voom pipeline
    fit = edgeR::voomLmFit(dge,design=data.frame(intercept = 1, label=1-(label<0)),)
    res = limma::eBayes(fit)
    
  } else if(limma_version == "limma_voom"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying RLE Normalization
    dge = edgeR::normLibSizes(dge,method="TMM")
    
    #Applying limma-voom pipeline
    fit = edgeR::voomLmFit(dge,design=data.frame(intercept = 1, label=1-(label<0)),)
    res = limma::eBayes(fit)
    
  } else if(limma_version == "limma_voom_RLE"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying RLE Normalization
    dge = edgeR::normLibSizes(dge,method="RLE")
    
    #Applying limma-voom pipeline
    fit = edgeR::voomLmFit(dge,design=data.frame(intercept = 1, label=1-(label<0)),)
    res = limma::eBayes(fit)
    
  } else if(limma_version == "limma_voom_TMMwsp"){
    
    #Creating DGE List Object
    dge = edgeR::DGEList(counts=df.test)
    
    #Applying TMMwsp Normalization
    dge = edgeR::normLibSizes(dge,method="TMMwsp")
    
    #Applying limma-voom pipeline
    fit = edgeR::voomLmFit(dge,design=data.frame(intercept = 1, label=1-(label<0)),)
    res = limma::eBayes(fit)
    
  } else {
    stop("limma version not recognised")
  }

  result=list(fit,res)
  names(result)=c("fit","res")
  
  return(result)
  
}