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


# select only relevant feature after RFE
train = train[,c("osevent","FUS","ATM","ERC1","TET1")]
test = test[,c("osevent","FUS","ATM","ERC1","TET1")]

# Cross validation
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

# Fitting SVM RBF to the Training set
svm_RBF <- train(osevent ~., data = train, method = "svmRadial",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)
# summary
svm_RBF

# test
test_pred <- predict(svm_RBF, newdata = test)
confusionMatrix(table(test_pred, test$osevent))
