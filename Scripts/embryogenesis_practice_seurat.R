#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(garnett)
library(SingleR)
library(celldex)
library(tidyverse)
library(pheatmap)



BiocManager::install("DESeq2")
BiocManager::install("SingleR")
BiocManager::install("celldex")

basepath <- "/Users/kaustubhgrama/Documents/GitHub/SCRI-KaustubhGrama/Scripts"

#load the data
cts <- ReadMtx(mtx = "/Users/kaustubhgrama/Desktop/Computer Science/R/Data/Small_Data/matrix.mtx.gz",
               features = "/Users/kaustubhgrama/Desktop/Computer Science/R/Data/Small_Data/features.tsv.gz",
               cells = "/Users/kaustubhgrama/Desktop/Computer Science/R/Data/Small_Data/barcodes.tsv.gz")

#create seurat object
org_data <- CreateSeuratObject(count = cts, project = "pmc3k", min.cells = 3, min.features = 200)

#QC-------------------------------------------------

#mitochondrial percentage
org_data[["percent.mt"]] <- PercentageFeatureSet(org_data, pattern = "^MT-")
#violin plot of the features, counts, and mitochondrial percentage
VlnPlot(org_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#subsetting the data
org_data <- subset(org_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalization---------------------------------------------------
org_data <- NormalizeData(org_data, normalization_method = "LogNormalize", scale.factor = 10000)

#Find Variable Features-------------------------------------------
org_data <- FindVariableFeatures(org_data, selection.method = "vst", nfeatures = 2000)
#isolating top10 variable features 
top10 <- head(VariableFeatures(org_data), 10)

#plot variable features with and without labels
plot1 <- VariableFeaturePlot(org_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#Scaling---------------------------------------------------
all.genes <- rownames(org_data)
org_data <- ScaleData(org_data, features = all.genes)

#PCA--------------------------------------------------------
org_data <- RunPCA(org_data, features = VariableFeatures(object = org_data))
DimHeatmap(org_data, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(org_data)

#Clustering--------------------------------
org_data <- FindNeighbors(org_data, dims = 1:15)
org_data <- FindClusters(org_data, resolution = 0.2)
org_data <- FindClusters(org_data, resolution = 0.5)
org_data <- FindClusters(org_data, resolution = 1.2)

#UMAP------------------------------------------------
org_data <- RunUMAP(org_data, dims = 1:15)

#visualization------------------------------------------

#view PCA clusters at different resolutions
DimPlot(org_data, reduction = "pca", group.by = "RNA_snn_res.0.2")
DimPlot(org_data, reduction = "pca", group.by = "RNA_snn_res.0.5")
DimPlot(org_data, reduction = "pca", group.by = "RNA_snn_res.1.2")

#view UMAP clusters at different resolutions
DimPlot(org_data, reduction = "umap", group.by = "RNA_snn_res.0.2")
DimPlot(org_data, reduction = "umap", group.by = "RNA_snn_res.0.5")
DimPlot(org_data, reduction = "umap", group.by = "RNA_snn_res.1.2")

#View UMAP by seurat clusters
seuratClusters.plot <- DimPlot(org_data, reduction = "umap", group.by = "seurat_clusters", label = TRUE)


#Find All Markers------------------------------------
all.markers <- FindAllMarkers(org_data, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = 'DESeq2', slot = 'counts')

#SINGLER-----------------------------------------------
#get reference data-------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

#run SingleR (default mode)
#default for SingleR is to perform annotation of each individual cell in the test dataset
org_data_counts <- GetAssayData(org_data, slot = 'counts')
main_pred = SingleR(test = org_data_counts, ref = ref, labels = ref$label.main)
fine_pred = SingleR(test = org_data_counts, ref = ref, labels = ref$label.fine)


#insert SingleR annotations into metadata
org_data$singleR.mainLabels <- main_pred$labels[match(rownames(org_data@meta.data), rownames(main_pred))]
org_data$singleR.fineLabels <- fine_pred$labels[match(rownames(org_data@meta.data), rownames(fine_pred))]



#View UMAP plot with SingleR annotations
singleR.mainLables.plot <- DimPlot(org_data, reduction = "umap", group.by = "singleR.mainLabels", label = TRUE)
singleR.fineLables.plot <- DimPlot(org_data, reduction = "umap", group.by = "singleR.fineLabels", label = TRUE)

