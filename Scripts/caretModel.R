#install packages
library(caret)
library(Seurat)
library(dplyr)
library(caretEnsemble)
install.packages("caretEnsemble")

seuobj89 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj89.RDS")
seuobj94 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj94.RDS")
seuobj110 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj110.RDS")
seuobj115 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj115.RDS")
seuobj125 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj125.RDS")

days <- c(seuobj89, seuobj94, seuobj110, seuobj115, seuobj125)

#preprocessing
for(day in days){
  print(day <- NormalizeData(day))
  print(day <- FindVariableFeatures(day))
  print(day <- ScaleData(day))
  print(day <- RunPCA(day))
}

#integrating all days into 1 object
achors <- FindIntegrationAnchors(object.list = days, dims = 1:30)
integrated_days <- IntegrateData(anchorset = achors, dims = 1:30)

integrated_days <- ScaleData(integrated_days)
integrated_days <- RunPCA(integrated_days)
integrated_days <- RunUMAP(integrated_days, dims = 1:30)


#Extract PCA embeddings
integrated_days_features <- Embeddings(integrated_days, "pca")
integrated_days_labels <- integrated_days@meta.data$Main_cluster_name

#create training and testing sets
set.seed(123)  # For reproducibility
trainIndex <- createDataPartition(integrated_days_labels, p = .8, list = FALSE, times = 1)
trainData <- integrated_days_features[trainIndex, ]
testData <- integrated_days_features[-trainIndex, ]
trainLabels <- integrated_days_labels[trainIndex]
testLabels <- integrated_days_labels[-trainIndex]

#train a model
# Define the training control
train_control <- trainControl(method = "cv", number = 10)
# Train the model (e.g., Random Forest)
model <- train(x = trainData, y = trainLabels, 
               method = "rf", 
               trControl = train_control)

models <- caretList(x = trainData, y = trainLabels, 
                    trControl = train_control,
                    methodList = c("rf", "svmRadial", "gbm"))
ensemble_model <- caretEnsemble(models)





#make predictions
test_predictions <- predict(model, testData)

testLabels <- factor(testLabels)
test_predictions <- factor(test_predictions, levels = levels(testLabels))


#confusion matrix
confusionMatrix(test_predictions, testLabels)

#make predictions on the integrated_day object
integrated_days_prediction <- predict(model, integrated_days_features)
integrated_days_prediction <- factor(integrated_days_prediction, levels = levels(testLabels))
integrated_days <- AddMetaData(object = integrated_days, metadata = integrated_days_prediction, col.name = "predicted_cell_type")



#make predictions on the seuobj115
seuobj115_features <- Embeddings(seuobj115, "pca")
seuobj115_prediction <- predict(model, seuobj115_features)
seuobj115 <- AddMetaData(object = seuobj115, metadata = seuobj115_prediction, col.name = "predicted_cell_type")
DimPlot(seuobj115, reduction = "umap", group.by = "predicted_cell_type", label = TRUE)








