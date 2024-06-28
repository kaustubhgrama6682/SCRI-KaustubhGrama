#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(garnett)


#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("cole-trapnell-lab/garnett")





basepath <- "/Users/kaustubhgrama/Desktop/Computer Science/R/"

#load the data
cts <- ReadMtx(mtx = "/Users/kaustubhgrama/Desktop/Computer Science/R/Small_Data/matrix.mtx.gz",
               features = "/Users/kaustubhgrama/Desktop/Computer Science/R/Small_Data/features.tsv.gz",
               cells = "/Users/kaustubhgrama/Desktop/Computer Science/R/Small_Data/barcodes.tsv.gz")

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

#view PCA clusters at different resolutions
DimPlot(org_data, reduction = "pca", group.by = "RNA_snn_res.0.2")
DimPlot(org_data, reduction = "pca", group.by = "RNA_snn_res.0.5")
DimPlot(org_data, reduction = "pca", group.by = "RNA_snn_res.1.2")

#view UMAP clusters at different resolutions
DimPlot(org_data, reduction = "umap", group.by = "RNA_snn_res.0.2")
DimPlot(org_data, reduction = "umap", group.by = "RNA_snn_res.0.5")
DimPlot(org_data, reduction = "umap", group.by = "RNA_snn_res.1.2")

#save seurat object
saveRDS(org_data, paste0(basepath, "org_data_seuratobj.RDS"))

#set object identity to chosen resolution
Idents(org_data) <- "RNA_snn_res.0.2"

#find all markers 
org_data.markers <- FindAllMarkers(org_data)

#Finding Differentially Expressed Features------------------------------------
org_data.markers %>% 
  group_by(cluster) %>% #grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% #filtering all groups by log2FC > 1
  slice_head(n = 10) %>% #selecting top 10 genes for each cluster
  ungroup() -> top10 #restaking each group (aka cluster ) into one variable dataframe

#heatmap of the top 10 upregulated genes per cluster res 0.2
DoHeatmap(org_data, features = top10$gene) + NoLegend() 

#Gene Ontology Analysis-------------------------------------------------------
cluster0_features <- top10$gene[1:10]
cluster1_features <- top10$gene[11:20]
cluster2_features <- top10$gene[21:30]
cluster3_features <- top10$gene[31:40]
cluster4_features <- top10$gene[41:50]
cluster5_features <- top10$gene[51:60]
cluster6_features <- top10$gene[61:70]
cluster7_features <- top10$gene[71:80]

FeaturePlot(org_data, features = cluster0_features, reduction = "umap")
FeaturePlot(org_data, features = cluster1_features, reduction = "umap")
FeaturePlot(org_data, features = cluster2_features, reduction = "umap")
FeaturePlot(org_data, features = cluster3_features, reduction = "umap")
FeaturePlot(org_data, features = cluster4_features, reduction = "umap")
FeaturePlot(org_data, features = cluster5_features, reduction = "umap")
FeaturePlot(org_data, features = cluster6_features, reduction = "umap")
FeaturePlot(org_data, features = cluster7_features, reduction = "umap")


#Create new seurat object with only cluster 0 cells
org_data.cluster0 <- subset(x = org_data, idents = "0")





