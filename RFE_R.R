set.seed(7)

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
mutations <- mutations[,-2] #remove id
mutations <- mutations[,-1] #remove id
data <- cbind(mutations, d1["osevent"])
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

# select only relevant feature after RFE
train = train[,c("osevent","FUS","ATM","ERC1","TET1")]
test = test[,c("osevent","FUS","ATM","ERC1","TET1")]

# train the model on selected features
rf = train(x=train[,-1], y=train[,1], method='rf')

# acc on train data
trainPred = predict(rf, train)
confusionMatrix(train[,1], trainPred)

testPred = predict(rf, test)
confusionMatrix(test[,1], testPred)
