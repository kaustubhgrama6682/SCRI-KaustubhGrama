# Load required libraries
library(readr)
library(dplyr)
library(Seurat)

folders <- c("GSM8012183_SIGAF2", "GSM8012186_SIGAG9", "GSM8012184_SIGAG11", "GSM8012185_SIGAG2", "GSM8012187_SIGAH10", "GSM8012188_SIGAH11", "GSM8012189_SIGAH2", "GSM8012190_SIGAH9", "GSM8012182_SIGAF11", "GSM8012181_SIGAF10", "GSM8012180_SIGAE1", "GSM8012179_SIGAD1", "GSM8012178_SIGAC1", "GSM8012177_SIGAB1")

seuobjs <- c()


for(folder in folders){
  basepath <- paste0("/Users/kaustubhgrama/Downloads/GSE252995_RAW/", folder, "/", folder, "_filtered_feature_bc_matrix_")

  #load the data
  cts <- ReadMtx(mtx = paste0(basepath, "matrix.mtx.gz"),
                 features = paste0(basepath, "features.tsv.gz"),
                 cells = paste0(basepath, "barcodes.tsv.gz"))
  
  seuobjs <- c(seuobjs,  CreateSeuratObject(count = cts, min.cells = 3, min.features = 200))
  
  
}

seuobjspreprocessed <- c()

for (seuobj in seuobjs){
  seuobj[["percent.mt"]] <- PercentageFeatureSet(seuobj, pattern = "^MT-")
  seuobj <- subset(seuobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  seuobj <- NormalizeData(seuobj, normalization.method = "LogNormalize", scale.factor = 10000)
  seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
  seuobj <- ScaleData(seuobj, features = rownames(seuobj))
  seuobj <- RunPCA(seuobj, features = VariableFeatures(object = seuobj))
  DimHeatmap(seuobj, dims = 1:18, cells = 500, balanced = TRUE)
  ElbowPlot(seuobj)
  seuobj <- FindNeighbors(seuobj, dims = 1:15)
  seuobj <- FindClusters(seuobj, resolution = 0.8)
  seuobj <- RunUMAP(seuobj, dims = 1:15)
  DimPlot(seuobj, reduction = "umap", group.by = "seurat_clusters")
  
  seuobjspreprocessed <- c(seuobjspreprocessed, seuobj)
  
}

seuobj1 <- seuobjspreprocessed[[1]]
subset <- seuobjspreprocessed[-c(1)]


merged_seuobj <- merge(x = seuobj1, y = subset)

merged_seuobj <- readRDS(file = "/Users/kaustubhgrama/Downloads/GSE252995_RAW/merged_seuobj")

merged_seuobj <- NormalizeData(merged_seuobj)
merged_seuobj <- FindVariableFeatures(merged_seuobj)
merged_seuobj <- ScaleData(merged_seuobj)
merged_seuobj <- RunPCA(merged_seuobj)
merged_seuobj <- FindNeighbors(merged_seuobj, dims = 1:15)
merged_seuobj <- FindClusters(merged_seuobj, resolution = 0.8)
merged_seuobj <- RunUMAP(merged_seuobj, dims = 1:15)
DimPlot(merged_seuobj, reduction = "umap")




