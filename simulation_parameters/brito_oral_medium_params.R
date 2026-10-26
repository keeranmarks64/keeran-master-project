##Writing Parameters to be sourced into generate_simulations file

#Filenames
baseline_meta_filename = "raw_data/meta_Brito_2016.tsv"
baseline_motu_filename = "raw_data/motus_Brito_2016.tsv"

#simulation parameters
ab.cutoff = 1e-04 #eliminates features that have less than this abundance (presumably across all samples)
prev.cutoff = 0.05 #eliminates features with less than this prevalence
log.n0 = 1e-05

# other parameters for the simulations
prop.markers = 0.1
class.balance = 0.5
subset_size = c(12,24,50,100,140)

# effect sizes
ab.scale = c(1, 1.25, 1.5, 2, 5, 10, 20)
prev.scale = c(0,0.1,0.2, 0.3)

#repeats
effect_comb_repeats = 15
subset_repeats = 25