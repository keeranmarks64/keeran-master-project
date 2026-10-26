import h5py
from pathlib import Path
import numpy as np
import pandas as pd
import sklearn as sk

# import sklearn.metrics
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.model_selection import LeaveOneOut, cross_val_predict, cross_val_score, RepeatedStratifiedKFold, cross_validate

def classifier_loocv(X, y, classifier):
    """
    Perform Leave-One-Out Cross-Validation using Random Forest Classifier.

    Parameters:
    X (array-like): Feature matrix of shape (n_samples, n_features).
    y (array-like): Target labels of shape (n_samples,).

    Returns:
    tuple: Mean accuracy, standard deviation of accuracy, and per-fold accuracies.
    """
    # Leave-one-out cross-validation to test the classifier
    loo = sk.model_selection.LeaveOneOut()

    if classifier == "random_forest_default":
        clf = sk.ensemble.RandomForestClassifier(
            n_estimators=100, #number of trees in the forest
            random_state=42, 
            n_jobs=-1, #number of jobs to run in parallel (-1 means using all available cores)
        )
    else:
        raise ValueError(f"Unsupported classifier: {classifier}")

    predict = sk.model_selection.cross_val_predict(clf, X, y, cv=loo, n_jobs=-1)

    return predict

def classifier_cv(X, y, classifier, n_splits, n_repeats):
    """
    Perform Cross-Validation using Random Forest Classifier.

    Parameters:
    X (array-like): Feature matrix of shape (n_samples, n_features).
    y (array-like): Target labels of shape (n_samples,).
    classifier (string): Classifier to use.
    n_splits (int): Number of splits for cross-validation.
    n_repeats (int): Number of repeats for cross-validation.

    """

    if classifier == "random_forest_default":
        clf = sk.ensemble.RandomForestClassifier(
            n_estimators=100, #number of trees in the forest
            random_state=42, 
            n_jobs=-1, #number of jobs to run in parallel (-1 means using all available cores)
        )
    else:
        raise ValueError(f"Unsupported classifier: {classifier}")


    cv = sk.model_selection.RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=42)

    cross_validation = sk.model_selection.cross_validate(clf, X, y, cv=cv, n_jobs=-1, scoring = ['roc_auc','precision','recall','accuracy'])

    return cross_validation



# Importing simulated data
simulation = "zeevi_classifier_test_tiny"
file_path = Path(f"datasets/{simulation}.h5")

f = h5py.File(file_path, "r")

#Listing Groups
groups = list(f.keys())
print(groups)

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
    print({"ab.scale": effect_combs.loc[group, "ab_scale"], "prev.scale": effect_combs.loc[group, "prev_scale"]})

#Preparing Table for storing the results of the classification metrics
results_table = pd.DataFrame(columns = ["Group","Group Name","Ab.scale","Prev.scale","Subset","Subset Repeat","Classifier","Mean Precision", "Mean Recall", " Mean AUROC", "Mean Accuracy"])

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

        precision_list = []
        recall_list = []
        auroc_list = []
        accuracy_list = []

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
            classifier = "random_forest_default"
            predict = classifier_loocv(X, y, classifier)

            true_labels = y
            predicted_labels = predict

            #computing accuracy metrics and appending them to list :
            precision_list.append(sklearn.metrics.precision_score(true_labels, predicted_labels))
            recall_list.append( sklearn.metrics.recall_score(true_labels, predicted_labels))
            auroc_list.append(sklearn.metrics.roc_auc_score(true_labels, predicted_labels))
            accuracy_list.append(sklearn.metrics.accuracy_score(true_labels, predicted_labels))

            results_table.loc[len(results_table.index)] = [group,effect_combs.loc[group, "group_name"],effect_combs.loc[group, "ab_scale"],effect_combs.loc[group, "prev_scale"],subset,subset_repeat,classifier, np.mean(precision_list), np.mean(recall_list), np.mean(auroc_list), np.mean(accuracy_list)]

filename = f"results_{simulation}.csv"
output_path = Path.cwd() / filename

results_table.to_csv(output_path, index=False)
