#load libraries
library(dplyr)
library(Seurat)
library(patchwork)

basepath <- "/Users/kaustubhgrama/Desktop/Computer Science/R/"


#load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/kaustubhgrama/Desktop/Computer Science/R/filtered_gene_bc_matrices/hg19")
#Initialize Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


# pbmc@assays$RNA
# View(pbmc@assays$RNA)

#QC----------------------------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Normalize Data---------------------------------------------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#identify 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling---------------------------------------------------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA-------------------------------------------------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1:5, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)

#Clustering----------------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:8)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- FindClusters(pbmc, resolution = 0.2)
pbmc <- FindClusters(pbmc, resolution = 1.2)

#UMAP-------------------------------------------------------------
pbmc <- RunUMAP(pbmc, dims = 1:8)

DimPlot(pbmc, reduction = "pca")
DimPlot(pbmc, reduction = "pca", group.by = "RNA_snn_res.0.5")
DimPlot(pbmc, reduction = "pca", group.by = "RNA_snn_res.0.2")
DimPlot(pbmc, reduction = "pca", group.by = "RNA_snn_res.1.2")

#save seurat object
saveRDS(pbmc, paste0(basepath, "pbmc_seuratobj.RDS"))

#set object identity to least resolution
Idents(pbmc) <- "RNA_snn_res.0.2"

#find all markers
pbmc.markers <- FindAllMarkers(pbmc)
View(pbmc.markers[pbmc.markers$cluster =="5",])

#Finding Differentially Expressed Features-------------------------
pbmc.markers %>%
  group_by(cluster) %>% #grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% #filtering all groups by log2FC > 1
  slice_head(n = 10) %>% #selecting top 10 genes for each cluster
  ungroup() -> top10 #restacking each group (aka cluster) into one variable dataframe i.e. "top10"

#heatmap of the top 10 upregulated genes per cluster res 0.2
DoHeatmap(pbmc, features = top10$gene) + NoLegend() 

#Gene Ontology Analysis------------------------------------------
naive_cd4_tcell_markers <- c("IL7R", "CCR7")
FeaturePlot(pbmc, features = naive_cd4_tcell_markers, reduction = "umap")

#Adding a module score based on the naive cd4 gene markers
pbmc <- AddModuleScore(pbmc, features = list("naive_cd4_tcell_markers"=naive_cd4_tcell_markers))

#renaming "cluster1" in the pbmc metadata
colnames(pbmc@meta.data)[9] <- "naive_cd4_tcell"

FeaturePlot(pbmc, features = "naive_cd4_tcell")




