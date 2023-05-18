# This script performs a classification (SVM, RF)
# on the joint dataset between mutations and
# clinical data (only on selected features for each dataset indipendently)

set.seed(7)

# load the library
library(mlbench)
library(caret)
library(e1071)
library(limma)
library("readxl")
library(dplyr)
library("dplyr")
library("faux")
library("DataExplorer")
library("caret")
library("randomForest")

data <- read.csv("./merged.csv")
# removing other labels and id's
data = subset(data, select = -c(X,AccessionNumber,dos,dpfs,X1stpfs.event) )
data$os.event = as.factor(data$os.event)

#split into train test
## 80% of the sample size
smp_size <- floor(0.80 * nrow(data))
train_ind <- sample(seq_len(nrow(data)), size = smp_size)
train <- data[train_ind, ]
test <- data[-train_ind, ]

# Cross validation
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

# Fitting SVM RBF to the Training set
svm_RBF <- train(os.event ~., data = train, method = "svmRadial",
                 trControl=trctrl,
                 preProcess = c("center", "scale"),
                 tuneLength = 10)
# summary
svm_RBF

# test
test_pred <- predict(svm_RBF, newdata = test)
confusionMatrix(table(test_pred, test$os.event))

#variable importance
geneImp=varImp(svm_RBF)
print(geneImp)
plot(geneImp)


# Random Forest
# train the model on selected features
# label is on col number 12
rf = train(x=train[,-12], y=train[,12], method='rf')

# acc on train data
trainPred = predict(rf, train)
confusionMatrix(train[,12], trainPred)

# test model
testPred = predict(rf, test)
confusionMatrix(test[,12], testPred)

#variable importance
geneImp=varImp(rf)
print(geneImp)
plot(geneImp)

