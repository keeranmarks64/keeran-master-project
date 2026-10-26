### Script to encapsulate apply_edgeR function

source("edgeR_versions.R")

apply_edgeR_group <- function(group,sim_implant,edgeR_version){
  #' Function to apply edgeR to a simbaR generated simulations, contained with the h5 file named
  #' 'sim_implant', specifically to the simulation specified in the string 'group'.
  #' 
  #' @param group: The specific name of the simulated group on which to apply edgeR, contained within the sim_implant file. 
  #' @param sim_implant: the name of the h5 file containing the simulated groups.
  #' @param edgeR_version: The specific version of the edgeR procedure to apply. 
  #' @return Returns a list containing p-values for each species, in each repetition, of each subset, of each simulated group.

  #extracting test data matrix from group
  complete_group = h5read(sim_implant, group)
  params <- h5read(file = sim_implant, name = "simulation_parameters")
  
  #Features matrix to analyse
  features = complete_group$features
  colnames(features) <- complete_group$sample_names
  rownames(features) <- complete_group$feature_names
  
  #Initialising list of pvalue matrices for each subset
  pvals_subsets = list()
  
  #Iterating through subsets
  subsets = names(complete_group$test_idx)
  for (subset in subsets) {
    
    #Getting Subset repeat indices
    subset_idx = complete_group[['test_idx']][[subset]]
    
    #Labelling each row for each repeat
    rownames(subset_idx) <- params$test_idx$names #labelling each repeat for each subset. 
    
    #Extracting group labels for subset repeats
    label=subset_idx["label",]
    
    #Subset repeat indices
    subset_repeats = subset_idx[which(grepl("^rep[0-9]*$", rownames(subset_idx))), 
                                , drop = FALSE]
    
    conf.mat=NULL #no confounders at this stage
    
    #Initialising p-value matrix 
    pvals=data.frame(row.names = complete_group$feature_names)
    
    #Iterating through repeats
    for(i in 1:dim(subset_repeats)[1]){
      
      #Dataset to be tested
      df.test = features[,subset_repeats[i,]] #taking ith repeat of subset
      
      #Normalising Data
      df.test = simbaR:::norm.data(df.test, norm="pass", log.n0) #normalisation step,  currently no normalisation
      
      #Applying edgeR, in different ways depending on version and extracting p-value matrix
      pvals_raw = apply_edgeR_version(df.test,label,edgeR_version)
      colname = paste0("rep_",i)
      
      pvals[[colname]] = pvals_raw
      
    }
    pvals_subsets[[subset]] <- as.matrix(pvals)
  }
  
  return(pvals_subsets)
}



