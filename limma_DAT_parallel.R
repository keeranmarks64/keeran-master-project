#Importing Simulated Data
library(simbaR)
library(purrr)
library(furrr)
library(tictoc)
source("limma_apply_parallel.R")

tic()
plan("multisession")
print(availableCores())
mb = 1000
options(future.globals.maxSize=mb*1024^2) 
args = commandArgs(trailingOnly = TRUE)

#Parameters for DAT
sim.implant_name = args[1] #simulation data filename
sim.implant_loc = paste0("datasets/",sim.implant_name,".h5")
sig = 0.05 #significance level of test
limma_versions = args[-1] #limma versions to a be applied

#Accessing group names
temp = h5ls(sim.implant_loc,recursive=FALSE)
groups = setdiff(temp$name, c('original_data', 'simulation_parameters')) #removing original data and simulation parameters from group names
sim_params = h5read(sim.implant_loc,"simulation_parameters")

## Differential Abundance Testing on Simulated Data 

for(limma_version in limma_versions){
  
  pvals_groups=future_map(groups,function(x) apply_limma_group(x,sim.implant_loc,limma_version)) 
  names(pvals_groups) = groups
  
  #Saving pvals to an h5 file
  h5save(pvals_groups,file=paste0('datasets/',sim.implant_name,"_",limma_version,'_pvals.h5'),name='pvals')
  
  results = future_map(groups, function(x){
    subset_results = list()
    
    for(subset in sim_params$test_idx$subsets){
      
      subset_name = paste0("subset_",as.character(subset))
      
      subset_results[[subset_name]]= eval.test( #BH correction applied here. 
        # full simulation file
        sim.location = sim.implant_loc,
        # which effect size combination
        group = x,
        # single p-value matrix
        res.mat = pvals_groups[[x]][[subset_name]],
        # desired significance threshold
        alpha = sig
        )
    }
    subset_results
    }) 
  names(results) = groups
  
  #saving results to file
  h5save(results,file=paste0('datasets/',sim.implant_name,"_",limma_version,'_results.h5'),name='results')
  
}

plan("sequential")

toc()

