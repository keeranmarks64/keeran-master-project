import sklearn as sk
import xgboost

def classifier_cv(X, y, classifier):
    """
    Perform Cross-Validation using Random Forest Classifier.

    Parameters:
    X (array-like): Feature matrix of shape (n_samples, n_features).
    y (array-like): Target labels of shape (n_samples,).
    classifier (string): Classifier to use.

    """

    if classifier == "random_forest_default":
        clf = sk.ensemble.RandomForestClassifier(
            n_estimators=100, #number of trees in the forest
            random_state=42, 
            n_jobs=-1 #number of jobs to run in parallel (-1 means using all available cores)
        )
    elif classifier == "decision_tree":
        clf = sk.tree.DecisionTreeClassifier(
            random_state=42
        )
    elif classifier == "logistic_regression":
        clf = sk.linear_model.LogisticRegression(
            random_state=42
        )
    elif classifier == "logistic_ridge_regression":
        clf = sk.linear_model.LogisticRegression(
            random_state=42, 
            l1_ratio=0
        )
    elif classifier == "logistic_lasso_regression":
        clf = sk.linear_model.LogisticRegression(
            random_state=42, 
            l1_ratio=1,
            solver="liblinear"
        )
    elif classifier == "svm_linear":
        clf = sk.svm.SVC(
            kernel="linear",
            random_state=42
        )
    elif classifier == "svm_rbf":
        clf = sk.svm.SVC(
            kernel="rbf",
            random_state=42 
        )
    elif classifier == "svm_poly":
        clf = sk.svm.SVC(
            kernel="poly",
            degree=3,
            random_state=42
        )
    elif classifier == "xgboost":
        clf = xgboost.XGBClassifier(
            random_state=42
        )


    return clf