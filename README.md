# keeran-master-project

### Data Collection Workflow for Simulated Microbiome Differential Abundance benchmarking
(Designed to be run from command line in a shell)
1.	Generating Simulated Dataset
Command: 
	R generate_simulations.R simulationName
simulationName is the name of the simulated dataset, which has a corresponding set of parameters contained with a file in the same directory named simulationName_params.R which is imported by the generate_simulations.R script. The script produces an h5 file containing the simulated data under the datasets folder called 
simulationName.h5
	The simulation parameters in the Wirbel study (Zeevi gut niche data) are
-	Abundance scales (ab.scale) : 1, 1.25, 1.5, 2, 5, 10, 20
-	Prevalence shift (prev.scale): 0,0.1,0.2, 0.3
-	Repetitions of each effect combination (effect_comb_repeats): 20
-	Sample Sizes (subset_size): 12, 24, 50, 100, 200, 400, 800
-	Sample Subset Repeats (subset_repeats): 50 
-	Abundance cutoff (ab.cutoff): 1e-04
-	Prevalence cutoff (prev.cutoff): 0.05
-	(log.n0): 1e-05
-	Proportion of Differentially Abundance species in the different groups (prop.markers): 0.1
-	Relative Proportion of samples (class.balance):
-	baseline_meta_filename: "raw_data/meta_Zeevi_2014.tsv"
-	baseline_motu_filename: "raw_data/motus_Zeevi_2014.tsv"

2.	Applying limma versions to simulated dataset and storing results
Command:
Rscript limma_DAT_parallel.R simulationName limmaVersion1 limmaVersion2 limmaVersion3
simulationName is the name of the simulated dataset that the limma versions will be applied to (no .h5 suffix). limmaVersionX are the names of the limma versions that will be applied to the dataset, with each producing its own set of p-values and results, stored under the datasets folder in files named:
	simulationName_limmaVersionX_pvals.h5
	simulationName_limmaVersionX_results.h5
3.	Applying EdgeR versions to simulated dataset and storing results:
	Command:
Rscript edgeR_DAT_parallel.R simulationName edgeRVersion1 edgeRVersion2 edgeRVersion3
simulationName is the name of the simulated dataset that the edgeR versions will be applied to (no .h5 suffix). edgeRVersionX are the names of the limma versions that will be applied to the dataset, with each producing its own set of p-values and results, stored under the datasets folder in files named:
	simulationName_ edgeRVersionX _pvals.h5
	simulationName_ edgeRVersionX _results.h5
4.	Plotting Results from files 
In the notebook plotting.Rmd you can run the functions plot_linegraphs_by_effect and generate_linegraphs_by_effect and generate_linegraphs_by_subset to reproduce the line graphs for plotting the accuracy metrics as a function of sample and effect size for a given simulated dataset and results. 


### Data Collection Workflow for Simulated Microbiome Classification Benchmarking

1. Generate Simulated Data as above: 

2. Apply classifier to simulated data with given classifiers with the following command:

  python3.13 classification_script.py simulationName classifier1 classifier2

3. Plot generated results using the classification_plotting.ipynb notebook and run the successive cells, changing the simulation name and set of classifiers to plot the relevant simulation/classifier plots.
