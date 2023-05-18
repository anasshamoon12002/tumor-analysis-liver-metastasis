# Tumor-Metastases-Analysis
Predict/Classify the survival (OS event) in colorectal cancer after liver surgery

In Classification folder we performed ML only on clinical data, and then on clinical data + mutation data (both reduced by the feature selection):
  1. SVC.ipynb and KNN_RF.ipynb classify only on clinical data (reduced with RFE) and they are in python
  2. RFE_R.R makes the recursive feature selection on mutations (mutations.csv) and select 4 important genes
  3. ClassificationD1_D2.R make classification on the merged dataset (merged.csv) -> reduced clinical data + reduced mutation data
