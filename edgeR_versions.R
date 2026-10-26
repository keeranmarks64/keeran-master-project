### Script to apply specific edgeR version to feature and label matrices

apply_edgeR_version <- function(df.test,label,edgeR_version){
  #' Function to apply a specific version of edgeR procedure to data and label matrices.
  #' 
  #' 
  #' @param df.test: data matrix of count values for each species and sample
  #' @param label: vector of length equivalent to the number of samples, representing the samples that are part of the case (1) and control (-1) groups. 
  #' @param edgeR_version: a string specifiying the specific version of the edgeR procedure to be applied. 
  #' @return returns the fitted edgeR object and p-value results matrix for the data provided. 
  
  #identifying only the positive count values and assigning their labels
  nonzero.idx = which(rowSums(df.test) > 0)
  f.counts = df.test[nonzero.idx, ]
  groups = rep("pos", ncol(f.counts))
  groups[label == -1] = "neg"
  names(groups) = colnames(df.test)
  
  #applying edgeR to non-zero count values
  y = edgeR::DGEList(counts = f.counts, group = groups)
  t <- exp(rowMeans(log(f.counts)))
  
  if(edgeR_version == "edgeR_v1"){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "TMM")
    }
    
    y = edgeR::estimateDisp(y)
    et = edgeR::exactTest(y)
    res = edgeR::topTags(et, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == "edgeR_v1_TMMwsp"){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "TMMwsp")
    }
    
    y = edgeR::estimateDisp(y)
    et = edgeR::exactTest(y)
    res = edgeR::topTags(et, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == "edgeR_v1_RLE"){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "RLE")
    }
    
    y = edgeR::estimateDisp(y)
    et = edgeR::exactTest(y)
    res = edgeR::topTags(et, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == 'edgeR_v2'){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "TMM")
    }    
    
    design=data.frame(intercept = -1, label)
    y = edgeR::estimateDisp(y,design=design)
    fit = edgeR::glmFit(y,design=design)
    lrt = edgeR::glmLRT(fit)
    res = edgeR::topTags(lrt, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == 'edgeR_v2_TMMwsp'){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "TMMwsp")
    }    
    
    design=data.frame(intercept = -1, label)
    y = edgeR::estimateDisp(y,design=design)
    fit = edgeR::glmFit(y,design=design)
    lrt = edgeR::glmLRT(fit)
    res = edgeR::topTags(lrt, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == 'edgeR_v2_RLE'){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "RLE")
    }    
    
    design=data.frame(intercept = -1, label)
    y = edgeR::estimateDisp(y,design=design)
    fit = edgeR::glmFit(y,design=design)
    lrt = edgeR::glmLRT(fit)
    res = edgeR::topTags(lrt, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == "edgeR_ql"){

    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "TMM")
    }   
    
    design=data.frame(intercept = -1, label)
    y = edgeR::estimateDisp(y,design=design)
    fit = edgeR::glmQLFit(y,design=as.matrix(design))
    qlf = edgeR::glmQLFTest(fit)
    res = edgeR::topTags(qlf, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == "edgeR_ql_TMMwsp"){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "TMMwsp")
    }   
    
    design=data.frame(intercept = -1, label)
    y = edgeR::estimateDisp(y,design=design)
    fit = edgeR::glmQLFit(y,design=as.matrix(design))
    qlf = edgeR::glmQLFTest(fit)
    res = edgeR::topTags(qlf, n = length(nonzero.idx), sort.by = "none")
    
  } else if(edgeR_version == "edgeR_ql_RLE"){
    
    if (any(t > 0)) {
      y = edgeR::normLibSizes(y, method = "RLE")
    }   
    
    design=data.frame(intercept = -1, label)
    y = edgeR::estimateDisp(y,design=design)
    fit = edgeR::glmQLFit(y,design=as.matrix(design))
    qlf = edgeR::glmQLFTest(fit)
    res = edgeR::topTags(qlf, n = length(nonzero.idx), sort.by = "none")
    
  } else {
    stop("edgeR version not recognised")
  }
  
  #Extracting p-values
  p.val = rep(1, nrow(df.test))
  names(p.val) <- rownames(df.test)
  p.val[nonzero.idx] = res$table[, "PValue"]
  
  return(p.val)
  
}
