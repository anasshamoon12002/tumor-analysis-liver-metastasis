import pandas as pd
import pickle
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, classification_report
from sklearn.linear_model import SGDClassifier
import numpy as np

np.random.seed(42)

# Load the dataset into a pandas DataFrame
df = pd.read_csv("dataset/mutationsComplete.csv", index_col=False)
df = df.drop("Unnamed: 0", axis=1)
# print(df.head())

# exit()

df = df.loc[:, ~df.columns.isin(['AccessionNumber', '1stpfs event', 'dpfs', 'dos'])]
# df = df[df["safety"]==1] #select safety analysis
# df = df.drop('PatientCode',axis=1) #drop one row with Nan value
# df = df[~df.isin([-99]).any(axis=1)] #drop any rows with -99 value
# # For now, remove all rows with Nan values
# # remove rows with NaN values, it's actually only one row
# df = df.dropna()

print(len(df))

train, test = train_test_split(df, test_size=0.2, random_state=42)

# Separate the target variable from the features
y_train = train['os event']
x_train = train.drop('os event', axis=1)
y_test = test['os event']
x_test = test.drop('os event', axis=1)

# RFE (Feature Selection)

print("\nRecursive Feature Selection\n")

min_features_to_select = 1
# RFE
rfecv = RFECV(
    estimator=RandomForestClassifier(),
    step=1, # number of feature to eliminate per iteration
    cv=StratifiedKFold(5),
    scoring="accuracy",
    min_features_to_select=min_features_to_select
)
# fit
rfecv.fit(x_train, y_train)

# print optimal number of feature selected
print(f"Optimal number of features: {rfecv.n_features_}")

print(x_train.columns[(rfecv.get_support())])

n_scores = len(rfecv.cv_results_["mean_test_score"])
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Mean test accuracy")
plt.errorbar(
    range(min_features_to_select, n_scores + min_features_to_select),
    rfecv.cv_results_["mean_test_score"],
    #yerr=rfecv.cv_results_["std_test_score"],
)
plt.title("Recursive Feature Elimination")
plt.show()

print(max(rfecv.cv_results_["mean_test_score"]))

print(rfecv.ranking_) # vote rank per feature

print(rfecv.support_) # selected features


# Random Forest prediction with selected data

print("\nRandom Forest on Selected Data\n")

# predict
y_pred = rfecv.predict(x_test)

# compute accuracy
print(accuracy_score(y_test, y_pred))


# Selected Data with Relevant Features

# data with only the selected features, converted automatically into narray
selected_x_train = rfecv.transform(x_train)
selected_x_test = rfecv.transform(x_test)

# KNN on selected data

print("\nK Nearest Neighbours on Selected Data\n")

# Creating the hyperparameter grid
param_grid = {'n_neighbors': [3,5,10,15,20],
              'weights': ['uniform','distance'],
              'metric': ['euclidean','minkowski','cosine']}

 
# Instantiating knn classifier
knn = KNeighborsClassifier()
 
# Instantiating the GridSearchCV object
knn_cv = GridSearchCV(knn, param_grid, cv = 5, refit=True)
 
knn_cv.fit(selected_x_train, y_train)
 
# Print the tuned parameters and score
print("Tuned KNN Parameters: {}".format(knn_cv.best_params_))
print("Best score on validation is {}".format(knn_cv.best_score_))

# Model Assessment

# prediction
grid_predictions = knn_cv.predict(selected_x_test)
# accuracy
print(accuracy_score(y_test, grid_predictions))
# print classification report 
print(classification_report(y_test, grid_predictions)) 