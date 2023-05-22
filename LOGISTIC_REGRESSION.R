library(stats)   # For base R functions
library(caret)

#load data
data <- read.csv("merged.csv")
#data = subset(data, select = -c(X,AccessionNumber, tstag, nstag, age..random, KDM5A, HIPPO, PIK3CG, PIK3R1, ATP1A1, FUS, WNT, TET1, RTK_RAS, TP53, CELL_CYCLE))
#data = subset(data, select = -c(X,AccessionNumber) )
data = subset(data, select = -c(X,AccessionNumber,dos,dpfs,X1stpfs.event) )


set.seed(42)
train_indices <- createDataPartition(data$os.event, p = 0.75, list = FALSE)
training_set <- data[train_indices, ]
test_set <- data[-train_indices, ]


model <- glm(os.event ~ ., data = training_set, family = binomial)
predictions <- predict(model, newdata = test_set, type = "response")

# Convert probabilities to class labels (using a threshold of 0.5)
predicted_labels <- ifelse(predictions > 0.5, 1, 0)

# Create a confusion matrix
true_labels <- test_set$os.event
confusion_matrix <- table(true_labels, predicted_labels)

# Calculate evaluation metrics
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
specificity <- confusion_matrix[1, 1] / sum(confusion_matrix[1, ])
f1_score <- 2 * precision * recall / (precision + recall)

# Calculate AUC-ROC
library(pROC)
roc_obj <- roc(true_labels, predictions)
auc_roc <- auc(roc_obj)

# Print the evaluation metrics
cat("Confusion Matrix:\n")
print(confusion_matrix)

cat("\nAccuracy:", accuracy)
cat("\nPrecision:", precision)
cat("\nRecall (Sensitivity):", recall)
cat("\nSpecificity:", specificity)
cat("\nF1-Score:", f1_score)
cat("\nAUC-ROC:", auc_roc)




feature_names <- colnames(training_set)

target_variable <- "os.event" 
if (target_variable %in% feature_names) {
  feature_names <- feature_names[feature_names != target_variable]
}

feature_names
coefficients <- abs(coef(model))
coefficients <- coefficients[-1]

# Create a data frame with feature names and corresponding coefficient magnitudes
feature_importance <- data.frame(Feature = feature_names, Importance = coefficients)

# Sort the data frame in descending order of importance
feature_importance <- feature_importance[order(feature_importance$Importance, decreasing = TRUE), ]

# Plot the feature importance
library(ggplot2)
ggplot(data = feature_importance, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Feature") +
  ylab("Importance") +
  ggtitle("Feature Importance in Logistic Regression") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
















# Assuming your logistic regression model is stored in the variable 'model'
# Assuming your training data is stored in the variable 'train_data'

# Predict the class labels on the training set
train_pred <- predict(model, training_set, type = "response")

# Convert the predicted probabilities to class labels (0 or 1)
train_pred_labels <- ifelse(train_pred > 0.5, 1, 0)

# Calculate the accuracy on the training set
train_accuracy <- mean(train_pred_labels == training_set$os.event)

# Print the accuracy
cat("Training Accuracy:", train_accuracy, "\n")

