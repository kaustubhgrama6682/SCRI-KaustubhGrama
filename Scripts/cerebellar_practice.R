#load libraries 
library(dplyr)
library(Seurat)
library(patchwork)

basepath <- "/Users/kaustubhgrama/Desktop/Computer Science/R/"

#load the dataset, initialize seurat object
cerebellar <- readRDS("/Users/kaustubhgrama/Desktop/Computer Science/R/cerebellar_embryogenesis_seuratobject.RDS")

#identify variable features
cerebellar <- FindVariableFeatures(cerebellar, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cerebellar), 10)

#plot variable features
plot1 <- VariableFeaturePlot(cerebellar)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling
all.genes <- rownames(cerebellar)
cerebellar <- ScaleData(cerebellar, features = all.genes)


#UMAP
cerebellar <- RunUMAP(cerebellar, dims = 1:15)
DimPlot(cerebellar, reduction = "umap")