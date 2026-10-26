#Importing Simulated Data
library(simbaR)
library(purrr)
library(furrr)
library(tictoc)
source("edgeR_apply_parallel.R")

tic()
plan("multisession")
#plan("sequential")
print(availableCores())
mb = 1000
options(future.globals.maxSize=mb*1024^2) 

args = commandArgs(trailingOnly = TRUE)

#Parameters for DAT
sim.implant_name = args[1] #simulation data filename
sim.implant_loc = paste0("datasets/",sim.implant_name,".h5")
sig = 0.05 #significance level of test
edgeR_versions = args[-1] #edgeR versions to a be applied

#Accessing group names
temp = h5ls(sim.implant_loc,recursive=FALSE)
groups = setdiff(temp$name, c('original_data', 'simulation_parameters')) #removing original data and simulation parameters from group names
sim_params = h5read(sim.implant_loc,"simulation_parameters")

## Differential Abundance Testing on Simulated Data 

for(edgeR_version in edgeR_versions){
  
  pvals_groups=future_map(groups,function(x) apply_edgeR_group(x,sim.implant_loc,edgeR_version),.progress=TRUE) 
  names(pvals_groups) = groups
  
  #Saving pvals to an h5 file
  h5save(pvals_groups,file=paste0('datasets/',sim.implant_name,"_",edgeR_version,'_pvals.h5'),name='pvals')
  
  ## Evaluating and Storing test results
  
  results = future_map(groups, function(x){
    subset_results = list()
    
    for(subset in sim_params$test_idx$subsets){
      
      subset_name = paste0("subset_",as.character(subset))
      
      subset_results[[subset_name]]= eval.test(
        # full simulation file
        sim.location = sim.implant_loc,
        # which effect size combination
        group = x,
        # single p-value matrix
        res.mat = pvals_groups[[x]][[subset_name]],
        # desired significance threshold
        alpha = sig)
    }
    subset_results
  },.progress=TRUE) 
  names(results) = groups
  
  #saving results to file
  h5save(results,file=paste0('datasets/',sim.implant_name,"_",edgeR_version,'_results.h5'),name='results')
  
}
toc()

