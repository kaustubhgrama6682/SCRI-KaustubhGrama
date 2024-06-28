#load libraries
library(dplyr)
library(Seurat)
library(patchwork)

basepath <- "/Users/kaustubhgrama/Desktop/Computer Science/R/"

#load the dataset
mouse <- readRDS("/Users/kaustubhgrama/Desktop/Computer Science/R/seurat_object_E9.5.rds")

#QC---------------------------------------------------------------------------
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^MT-")
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)s
mouse <- subset(mouse, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize----------------------------------------------------------------------
mouse <- NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)

#Identify variable features-----------------------------------------------------
top10 <- head(VariableFeatures(mouse), 10)

#plot variable features
plot1 <- VariableFeaturePlot(mouse)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling--------------------------------------------------------------------
all.genes <- rownames(mouse)
mouse <- SCTransform(mouse, verbose = FALSE)