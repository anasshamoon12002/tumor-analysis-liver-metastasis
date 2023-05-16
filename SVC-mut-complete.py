import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.feature_selection import SelectKBest, chi2, f_classif
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, cross_val_score, StratifiedShuffleSplit
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import pandas as pd
import pickle
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.feature_selection import RFE
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import SGDClassifier

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

min_features_to_select = 1
# RFE
rfe = RFECV(
    estimator=RandomForestClassifier(),
    step=1, # number of feature to eliminate per iteration
    cv=StratifiedKFold(5),
    scoring="accuracy",
    min_features_to_select=min_features_to_select
)
# fit
rfe.fit(x_train, y_train)

# print optimal number of feature selected
print(f"Optimal number of features: {rfe.n_features_}")

print(x_train.columns[(rfe.get_support())])

print(max(rfe.cv_results_["mean_test_score"]))

# Get the selected features
x_train_selected = rfe.transform(x_train)
x_test_selected = rfe.transform(x_test)

# Define the parameter grid for GridSearchCV
param_grid = {
    'C': [0.8, 0.9, 1, 2, 3, 4, 5, 6], 
    'kernel': ['rbf', 'sigmoid'], 
    'gamma': [0.000001, 0.000002, 0.000003, 0.000004, 0.000005]
}

# Perform hyperparameter tuning with cross-validation
cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
grid_search = GridSearchCV(SVC(random_state=42), param_grid, cv=cv, scoring='accuracy')
grid_search.fit(x_train_selected, y_train)

# Print the best hyperparameters
print("Best hyperparameters:", grid_search.best_params_)

# Evaluate the performance of the best model using cross-validation
cv_results = cross_val_score(grid_search.best_estimator_, x_train_selected, y_train, cv=cv, scoring='accuracy')
print("Mean accuracy on val:", np.mean(cv_results))

# Train the final model using the selected features and best hyperparameters
clf = SVC(**grid_search.best_params_, random_state=42)
clf.fit(x_train_selected, y_train)

# Predict the labels of the dataset using the final model
y_pred = clf.predict(x_test_selected)

# Calculate the accuracy of the final model
accuracy = accuracy_score(y_test, y_pred)
print('Accuracy on test:', accuracy)