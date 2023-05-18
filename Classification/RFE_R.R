# this script perform feature selection on mutation data,
# performing more runs on different seed and then see which genes
# appears more in the list

set.seed(876)

# load the library
library(mlbench)
library(caret)
library(limma)
library("readxl")
library(dplyr)
library("dplyr")
library("faux")
library("DataExplorer")
library("caret")
library("randomForest")

# read the two datasets
d1 <- read_excel("./d1.xlsx")
mutations <- read.csv("./mutations.csv")
# remove id's
mutations <- mutations[,-2] 
mutations <- mutations[,-1] 
# bind label and mutations dataset
data <- cbind(mutations, d1["osevent"])
# convert into factor
data$osevent = as.factor(data$osevent)

#split into train test
## 80% of the sample size
smp_size <- floor(0.80 * nrow(data))
train_ind <- sample(seq_len(nrow(data)), size = smp_size)
train <- data[train_ind, ]
test <- data[-train_ind, ]

# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(train[,-569], train[,569], metric="Accuracy", rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))


# ------Can exhecute starting on this line------

# select only relevant feature after RFE (IN COMMON: did it already in python)
# "FUS","ATM","ERCC3","ERC1" it's the result of the python script which
# check the ranking of each word in a string list
train = train[,c("osevent","FUS","ATM","ERCC3","ERC1")]
test = test[,c("osevent","FUS","ATM","ERCC3","ERC1")]

# train RF model on selected features
rf = train(x=train[,-1], y=train[,1], method='rf')

# acc on train data
trainPred = predict(rf, train)
confusionMatrix(train[,1], trainPred)

# test model
testPred = predict(rf, test)
confusionMatrix(test[,1], testPred)

#variable importance
geneImp=varImp(rf)
print(geneImp)
plot(geneImp)
