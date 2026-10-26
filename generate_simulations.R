## Import + Setup
#Importing packages for simulation
library(tidyverse) # for data wrangling
library(simbaR)    #SimbaR for signal implantation

args = commandArgs(trailingOnly = TRUE) #only argument is name of simulation to be set
source(paste0("simulation_parameters/",args[1],"_params.R")) #getting simulation parameters 

#Importing Data
meta = read_tsv(baseline_meta_filename, col_types = cols()) %>%
  filter(!is.na(Sample_ID)) %>%
  as.data.frame()
rownames(meta) <- meta$Sample_ID

feat=read.table(baseline_motu_filename, sep='\t', stringsAsFactors = FALSE,
                      check.names = FALSE, quote = '', comment.char = '',
                      row.names = 1, header = TRUE)

##Signal Implantation
# give the simulation a filename
sim.implant <- paste0('datasets/',args[1],'.h5')

# create simulations
create.data.simulation(feat = feat,
                       meta = meta,
                       sim.location = sim.implant,
                       sim.type ='cross-section',
                       sim.method ='resampling',
                       filt.params = list(ab.cutoff=as.numeric(ab.cutoff),
                                          prev.cutoff=prev.cutoff,
                                          log.n0=as.numeric(log.n0)),
                       sim.params = list(ab.scale=ab.scale,
                                         prev.scale=prev.scale,
                                         prop.markers=prop.markers,
                                         conf=NULL,
                                         class.balance=class.balance,
                                         feature.type='all',
                                         repeats=effect_comb_repeats))

# create testing indices
create.test.idx(sim.location=sim.implant,
                subsets=subset_size,
                repetitions=subset_repeats)
