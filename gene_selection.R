data=read.csv("C:\\Users\\botru\\Desktop\\mutationsComplete.csv")
library("readxl")
library(caret)
library(glmnet)
data2=read_excel("C:\\Users\\botru\\Desktop\\d1.xlsx")
data2<-data2[7:69]
data2
print(colnames(data2))
library(dplyr)
data2 <- select_if(data2, is.numeric)
data2 <- data2[, !(names(data2) %in% c("oors", "dpfs", "1stpfs event","dos"))]


data <- data[,-1]
data <- data[,-3]
data <- data[,-3]
data <- data[,-3]
colnames(data)
ncol(data)

####################### PARTITION TEST SET & TRAINING SET #####################

library(caret)

# Imposta un seed per la riproducibilità
set.seed(42)

# Dividi il dataset in training set e test set
train_indices <- createDataPartition(data$os.event, p = 0.75, list = FALSE)
training_set <- data[train_indices, ]
test_set <- data[-train_indices, ]

training_patients = training_set[["AccessionNumber"]]
test_patients = test_set[["AccessionNumber"]]

training_set <- training_set[, -1]
test_set <- test_set[, -1]

nrow(training_set)
nrow(test_set)

###############################################################################

####################### DIFFERENTIAL EXPRESSION ###################################


library(GEOquery)
library(limma)

alive = training_set[["os.event"]]==0
alive<-alive*1
alive

dead = training_set[["os.event"]]==1
dead<-dead*1
dead

design2=cbind(ALIVE=alive,DEAD=dead)
design2

fitMy=lmFit(t(training_set),design2)
cont.matrix <- makeContrasts(DEAD-ALIVE, levels=design2)

fit2 <- contrasts.fit(fitMy, cont.matrix)
fit2 <- eBayes(fit2)

de=topTable(fit2,adjust='BH',number=100)
plot(de$logFC,-log(de$P.Value,10),xlim=c(-5,5))


######################## EXTRACT MOST IMPORTANT GENES ######################

# Ordina la tabella dei risultati per i valori di p corretti per il multiplo confronto
de_sorted <- de[order(de$P.Value, decreasing = FALSE), ]
de_sorted


# Specifica il valore di soglia per il P-value e il fold change
pvalue_threshold <- 0.05
fold_change_threshold <- 1.5

# Seleziona i geni con un P-value inferiore alla soglia e un fold change superiore alla soglia
significant_genes <- de_sorted[de_sorted$P.Value < pvalue_threshold & abs(de_sorted$logFC) > log2(fold_change_threshold), ]

# Ordina la tabella dei risultati filtrata per il P-value
significant_genes <- significant_genes[order(significant_genes$P.Value), ]
significant_genes
important_genes <- rownames(significant_genes)
important_genes <- important_genes[-1]
length(important_genes)


important_genes <- important_genes[1:20]
important_genes

###################################################################













########################## LOGISTIC REGRESSION ###############################
library(stats)

selected_genes <- training_set[, important_genes]
train_data <- cbind(selected_genes, os.event = training_set$os.event)
model <- glm(os.event ~ ., data = train_data, family = binomial)
test_data <- test_set[, important_genes]
predictions <- predict(model, newdata = test_data, type = "response")

# Trasforma le probabilità in classi (0 o 1) utilizzando una soglia di decisione (ad esempio, 0.5)
predicted_classes <- ifelse(predictions >= 0.5, 1, 0)

# Calcola l'accuracy
accuracy <- sum(predicted_classes == test_set$os.event) / length(test_set$os.event)

# Stampa l'accuracy
print(paste("Accuracy:", accuracy))

#0.7033

###################################### MERGE CLINICAL DATA AND GENES #####



all =read.csv("C:\\Users\\botru\\Desktop\\all.csv")
print(colnames(all))
all <- all[, !(names(all) %in% c("X", "AccessionNumber", "os.event_y", "Unnamed..0","dos","dpfs","X1stpfs.event"))]
library(dplyr)
all <- all[1:61]

# Rimuovi le colonne non numeriche dal dataframe "df"
all <- all %>% select_if(is.numeric)



# Imposta un seed per la riproducibilità
set.seed(42)

# Dividi il dataset in training set e test set
train_indices2 <- createDataPartition(all$"os.event_x", p = 0.75, list = FALSE)
training_set2 <- all[train_indices2, ]
test_set2 <- all[-train_indices2, ]
training_patients2 = training_set2[["AccessionNumber"]]
test_patients2 = test_set2[["AccessionNumber"]]
training_set2 <- training_set2[, !(names(training_set2) %in% c("AccessionNumber"))]
test_set2 <- test_set2[, !(names(test_set2) %in% c("AccessionNumber"))]

nrow(training_set2)
nrow(test_set2)



#install.packages("glmnet")


# Converti i tuoi dati nel formato richiesto da glmnet
print(colnames(training_set2))
training_set2
x <- as.matrix(training_set2[, !(names(training_set2) %in% c("os.event_x"))])
y <- as.numeric(training_set2$"os.event_x")  # Variabile di risposta

lasso_model <- glmnet(x, y, family = "binomial", alpha = 0)
lambda_min <- lasso_model$lambda[which.min(lasso_model$dev)]
lasso_coef <- coef(lasso_model, s = lambda_min)

# Estrai i nomi delle feature dalla matrice dei coefficienti
feature_names <- dimnames(lasso_coef)[[1]]

# Estrai i valori dei coefficienti
coefficient_values <- as.numeric(lasso_coef)

# Crea un data frame con i nomi delle feature e i valori dei coefficienti
coef_df <- data.frame(feature = feature_names, coefficient = coefficient_values)

# Rimuovi le righe con valore "."
coef_df <- coef_df[coef_df$coefficient != ".", ]

# Ordina il data frame in base ai valori dei coefficienti (in modo decrescente)
coef_df <- coef_df[order(-coef_df$coefficient), ]

# Estrai i nomi delle feature più importanti
important_features <- coef_df$feature
important_features
# Seleziona le prime 10 feature più importanti
top_features <- head(important_features, 10)

print(top_features)





# Effettua delle predizioni sul set di test
x_test <- as.matrix(test_set2[, !(names(test_set2) %in% c("os.event_x"))])
y_test <- test_set2$"os.event_x"
predictions <- predict(lasso_model, newx = x_test, type = "response")

# Valuta le prestazioni del modello
threshold <- 0.5
predicted_classes <- ifelse(predictions > threshold, 1, 0)
accuracy <- mean(predicted_classes == y_test)
print(accuracy)

########################### USING ONLY TOP FEATURES ############################


library(stats)

top_features<-top_features[-1]

selected_genes <- training_set2[, top_features]
train_data <- cbind(selected_genes, os.event_x = training_set2$os.event_x)
model <- glm(os.event_x ~ ., data = train_data, family = binomial)
test_data <- test_set2[, top_features]
predictions <- predict(model, newdata = test_data, type = "response")

# Trasforma le probabilità in classi (0 o 1) utilizzando una soglia di decisione (ad esempio, 0.5)
predicted_classes <- ifelse(predictions >= 0.5, 1, 0)
# Calcola l'accuracy
accuracy <- sum(predicted_classes == test_set2$os.event_x) / length(test_set2$os.event_x)

# Stampa l'accuracy
print(paste("Accuracy:", accuracy))






