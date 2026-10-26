import h5py
from pathlib import Path
import numpy as np
import pandas as pd
import sklearn as sk
import sys
from classifiers import classifier_cv

# Parsing command line arguments
simulation = sys.argv[1]
classifiers = sys.argv[2:]

print(f"Simulation: {simulation}")
print(f"Classifier(s): {classifiers}")

# Importing simulated data
file_path = Path(f"datasets/{simulation}.h5")
f = h5py.File(file_path, "r")

#Listing Groups
groups = list(f.keys())

#Getting effect combinations and associating them to group names
ab_scale= list(f["simulation_parameters"]['sim_params']['ab.scale'])
prev_scale = list(f["simulation_parameters"]['sim_params']['prev.scale'])
no_group_repeats = list(f["simulation_parameters"]['sim_params']['repeats'])

ab_scale_groups = pd.Series(np.repeat(ab_scale, len(prev_scale))).repeat(no_group_repeats)
prev_scale_groups = pd.Series(np.tile(prev_scale, len(ab_scale))).repeat(no_group_repeats)

groups = groups[0:len(groups)-2]
effect_combs = pd.DataFrame({"ab_scale": ab_scale_groups, "prev_scale": prev_scale_groups, "group_name": np.zeros(len(ab_scale_groups)).astype(str),"group":groups})
effect_combs = effect_combs.set_index("group")

for group in groups:
    effect_combs.loc[group,"group_name"] =  f"ab{effect_combs.loc[group, "ab_scale"]}prev{effect_combs.loc[group, "prev_scale"]}"

for classifier in classifiers:
    #Preparing Table for storing the results of the classification metrics
    results_table = pd.DataFrame(columns = ["Group","Group Name","Ab.scale","Prev.scale","Subset","Subset Repeat","Classifier","Precision", "Recall", "AUROC", "Accuracy", "Fit Time"])

    #Iterating through the groups and aggregating results from classifiers 
    for group in groups:
        print(group)
        group_keys = list(f[group].keys())

        #Extracting Feature Matrix as a dataframe
        features = pd.DataFrame(f[group][group_keys[1]], #Count Matrix 
                                index = f[group][group_keys[4]], #Sample IDs
                                columns = f[group][group_keys[0]] #Feature IDs
                                ).transpose() 

        #Extracting the group labels and DA species 
        labels = np.array(f[group][group_keys[2]])
        marker_idx = np.array(f[group][group_keys[3]])

        #Extracting Subset Indices
        subsets = list(f[group][group_keys[5]])

        #Iterate through the subsets
        for subset in subsets:
            subset_idx = np.array(f[group][group_keys[5]][subset])
            print(subset)

            #Iterate through the subset repeats
            for i in range(1, np.shape(subset_idx)[1]):
                subset_repeat = str(i)
                print("Repeat: ", subset_repeat)
                subset_repeat_idx = subset_idx[:,i]-1 #indices of subsets

                #Extracting the subset features and labels
                subset_repeat_labels = subset_idx[:,0]
                subset_features = features.iloc[:, subset_repeat_idx]

                # Prepare the data for classification
                X = subset_features.T #Sample feature matrix of dimension (n_samples, n_features)
                y = np.asarray(subset_repeat_labels).astype(int) #Sample labels of dimension (n_samples,)
                if classifier=='xgboost':
                    y = np.where(y == -1, 0, 1)

                clf = classifier_cv(X, y, classifier)

                cv = sk.model_selection.RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)

                cross_validation = sk.model_selection.cross_validate(clf, X, y, cv=cv, n_jobs=-1, scoring = ['roc_auc','precision','recall','accuracy'],)

                #computing accuracy metrics and appending them to list :
                precision = cross_validation['test_precision'].mean()
                recall = cross_validation['test_recall'].mean()
                auroc = cross_validation['test_roc_auc'].mean()
                accuracy = cross_validation['test_accuracy'].mean()
                fit_time = cross_validation['fit_time'].mean()

                results_table.loc[len(results_table.index)] = [group,effect_combs.loc[group, "group_name"],effect_combs.loc[group, "ab_scale"],effect_combs.loc[group, "prev_scale"],subset,subset_repeat,classifier, precision, recall, auroc, accuracy, fit_time]

    filename = f"classification_results/results_{simulation}_{classifier}.csv"
    output_path = Path.cwd() / filename

    results_table.to_csv(output_path, index=False,)
