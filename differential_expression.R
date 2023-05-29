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

training_set

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
de
#plot(de$logFC,de$P.Value,xlim=c(-5,5), pch=labels, cex=1)

# Plot genes based on log-fold change and p-value
plot(de$logFC, -log10(de$P.Value), xlim = c(-5, 4), ylim = c(0, 3), pch = 16, cex = 0.8,
     xlab = "Log-fold change", ylab = "-log10(p-value)")

# Add significance threshold lines
abline(h = -log10(0.05), col = "red", lty = 2)
abline(v = c(-1, 1), col = "blue", lty = 2)

# Add gene labels
text(de$logFC, -log10(de$P.Value), labels = rownames(de), pos = 4, cex = 0.8)



######################## EXTRACT MOST IMPORTANT GENES ######################

# Sort by pValue
de_sorted <- de[order(de$P.Value, decreasing = FALSE), ]
de_sorted


# Set thresholds
pvalue_threshold <- 0.05
fold_change_threshold <- 1

# Select genes with P-value under threshold and fold change over the threshold
significant_genes <- de_sorted[de_sorted$P.Value < pvalue_threshold & abs(de_sorted$logFC) > fold_change_threshold, ]

# Order by PValue
significant_genes <- significant_genes[order(significant_genes$P.Value), ]
significant_genes
important_genes <- rownames(significant_genes)
important_genes <- important_genes[-1]
length(important_genes)


important_genes <- important_genes[1:50]
important_genes

###################################################################








