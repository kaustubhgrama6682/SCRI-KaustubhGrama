#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(sctransform)

seuobj_atrt@active.assay = "SCT"
seuobj_dipg@active.assay = "SCT"
seuobj_mb@active.assay = "SCT"
seuobj89@active.assay = "SCT"
seuobj94@active.assay = "SCT"
seuobj110@active.assay = "SCT"
seuobj115@active.assay = "SCT"
seuobj125@active.assay = "SCT"

#developmental days
seuobj89 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj89.RDS")
seuobj94 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj94.RDS")
seuobj110 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj110.RDS")
seuobj115 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj115.RDS")
seuobj125 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj125.RDS")

#integration of all days

days <- c(
          "115" = seuobj115,
          "125" = seuobj125)

features <- SelectIntegrationFeatures(object.list = days)
anchors <- FindIntegrationAnchors(object.list = days, anchor.features = features)


seuobj_115_125 <- IntegrateData(anchorset = anchors)

seuobj_115_125 <- NormalizeData(seuobj_115_125)
seuobj_115_125 <- FindVariableFeatures(seuobj_115_125)

seuobj_115_125 <- ScaleData(seuobj_115_125)
seuobj_115_125 <- RunPCA(seuobj_115_125)
seuobj_115_125 <- RunUMAP(seuobj_115_125, dims = 1:30)
seuobj_115_125 <- FindNeighbors(seuobj_115_125, dims = 1:30)
seuobj_115_125 <- FindClusters(seuobj_115_125)

saveRDS(seuobj_115_125, file = "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/seuobj_115_125.RDS")


#tumors
seuobj_atrt <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj_atrt.RDS")
seuobj_dipg <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj_dipg.RDS")
seuobj_mb <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj_mb.RDS")


#atrt pre processing
seuobj_atrt[["percent.mt"]] <- PercentageFeatureSet(seuobj_atrt, pattern = "^MT-")
seuobj_atrt <- subset(seuobj_atrt, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seuobj_atrt <- NormalizeData(seuobj_atrt, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_atrt <- FindVariableFeatures(seuobj_atrt, selection.method = "vst", nfeatures = 2000)
seuobj_atrt <- ScaleData(seuobj_atrt, features = rownames(seuobj_atrt))
seuobj_atrt <- RunPCA(seuobj_atrt, features = VariableFeatures(object = seuobj_atrt))
DimHeatmap(seuobj_atrt, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj_atrt)
seuobj_atrt <- FindNeighbors(seuobj_atrt, dims = 1:15)
seuobj_atrt <- FindClusters(seuobj_atrt, resolution = 0.8)
seuobj_atrt <- RunUMAP(seuobj_atrt, dims = 1:15)
DimPlot(seuobj_atrt, reduction = "umap", group.by = "seurat_clusters")


#dipg pre processing
seuobj_dipg[["percent.mt"]] <- PercentageFeatureSet(seuobj_dipg, pattern = "^MT-")
seuobj_dipg <- subset(seuobj_dipg, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seuobj_dipg <- NormalizeData(seuobj_dipg, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_dipg <- FindVariableFeatures(seuobj_dipg, selection.method = "vst", nfeatures = 2000)
seuobj_dipg <- ScaleData(seuobj_dipg, features = rownames(seuobj_dipg))
seuobj_dipg <- RunPCA(seuobj_dipg, features = VariableFeatures(object = seuobj_dipg))
DimHeatmap(seuobj_dipg, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj_dipg)
seuobj_dipg <- FindNeighbors(seuobj_dipg, dims = 1:15)
seuobj_dipg <- FindClusters(seuobj_dipg, resolution = 0.8)
seuobj_dipg <- RunUMAP(seuobj_dipg, dims = 1:15)
DimPlot(seuobj_dipg, reduction = "umap", group.by = "seurat_clusters")


#SCTransform
# run sctransform
seuobj_atrt <- SCTransform(seuobj_atrt, vars.to.regress = "percent.mt", verbose = FALSE)
seuobj_dipg <- SCTransform(seuobj_dipg, vars.to.regress = "percent.mt", verbose = FALSE)


#this didn't work
####seuobj.list <- SCTransform(seuobj.list, vars.to.regress = "percent.mt", verbose = FALSE)

seuobj_int <- NormalizeData(seuobj_int, normalization.method = "LogNormalize", scale.factor = 10000)
seuobj_int <- FindVariableFeatures(seuobj_int, selection.method = "vst", nfeatures = 2000)
seuobj_int <- RunFastMNN(object.list = SplitObject(seuobj_int, split.by = "orig.ident"))
seuobj_int <- ScaleData(seuobj_int, features = rownames(seuobj_int))

seuobj_int <- RunPCA(seuobj_int, features = VariableFeatures(object = seuobj_int))
seuobj_int <- FindNeighbors(seuobj_int, dims = 1:15)
seuobj_int <- FindClusters(seuobj_int, resolution = 0.8)
seuobj_int <- RunUMAP(seuobj_int, dims = 1:15)






#Integration of tumors===========
#creating tumor seuobj list
seuobj.list <- list("atrt" = seuobj_atrt,
                    "dipg" = seuobj_dipg,
                    "mb" = seuobj_mb,
                    "s89" =  seuobj89, 
                    "s94" = seuobj94, 
                    "s115" = seuobj115,
                    "s110" = seuobj110, 
                    "s125" = seuobj125)

#identifying features
features <- SelectIntegrationFeatures(object.list = seuobj.list)
#creating anchors
anchors <- FindIntegrationAnchors(object.list = seuobj.list, anchor.features = features)
#creating integrated obj
seuobj_int <- IntegrateData(anchorset = anchors)



#switch to RNA
for(obj in seuobj.list){
  print(obj@active.assay <- "SCT")

}



#FeaturePlots for topics

topics <- c("Day125_Topic_33", "Day115_Topic_12", "Day110_Topic_14", "Day94_Topic_27",  "Day89_Topic_14")

for(topic in topics){
  print(FeaturePlot(seuobj_int, features = topic, reduction = "umap"))
}


#DEG
seuobj_int_markers <- FindAllMarkers(seuobj_int)
seuobj_int_markers %>%
  group_by(cluster) %>% #grouping dataframe by cluster column
  dplyr::filter(avg_log2FC > 1) %>% #filtering all groups by log2FC > 1
  slice_head(n = 10) %>% #selecting top 10 genes for each cluster
  ungroup() -> top10_int #restacking each group (aka cluster) into one variable dataframe i.e. "top10"

cluster1_genes <- top10_int$gene[11:20]
cluster2_genes <- top10_int$gene[21:30]
cluster0_genes <- top10_int$gene[1:10]
cluster8_genes <- top10_int$gene[81:90]
cluster17_genes <- top10_int$gene[171:180]
cluster20_genes <- top10_int$gene[201:210]
cluster18_genes <- top10_int$gene[181:190]
cluster23_genes <- top10_int$gene[231:240]



seuobj1_2_8 <- subset(seuobj_int, subset = seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 8)

seuobj0_17_18_20_23 <- subset(seuobj_int, subset = seurat_clusters == 0 | seurat_clusters == 17 | seurat_clusters == 18 | seurat_clusters == 20 | seurat_clusters == 23)

saveRDS(seuobj0_17_18_20_23, file = "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/seuobj0_17_18_20_23.RDS")

Day89_Topic_14_top15 <- c("NLGN1", "LSAMP", "ADAMTS6", "NFIB", "NFIA", "LINC00669", "RP11-649A16.1", "SYNE2", "EGFEM1P", "DIAPH3", "DLEU2", "WWOX", "CENPF", "ZBTB20", "RAD51B")
Day94_Topic_27_top15 <- c("SYNE2", "LINC00669", "EGFEM1P", "DIAPH3", "CENPF", "DLEU2", "MKI67", "FBXL7", "NFIA", "CENPP", "NLGN1", "TOP2A", "CEP128", "RFC3", "ASPM")
Day110_Topic14_top15 <- c("LINC00669", "SYNE2", "EGFEM1P", "NLGN1", "DIAPH3", "CENPP", "DLEU2", "WWOX", "RFC3", "APOLD1", "TOP2A", "CENPF", "BRIP1", "FBXL7", "MKI67")
Day115_Topic_12 <- c("LRRTM4", "PTN", "DACH1", "SLC1A3", "TNC", "PTPRM", "GPM6B", "RP11-649A16.1", "SEMA6A", "CTNNA3", "EGFEM1P", "CDH4", "NCKAP5", "NAALADL2", "DGKB")
Day125_Topic_33 <- c("LINC00669", "LSAMP", "NLGN1", "DCC", "DIAPH3", "CENPP", "DLEU2", "RFC3", "MKI67", "APOLD1", "WWOX", "TOP2A", "CENPF", "RAD51B", "EGFEM1P")

top15s <- c(Day89_Topic_14_top15, Day94_Topic_27_top15, Day110_Topic14_top15, Day115_Topic_12, Day125_Topic_33)
top15names <- c("Day89_Topic_14_top151", "Day94_Topic_27_top151","Day110_Topic14_top151", "Day115_Topic_121", "Day125_Topic_331")
top15namesnew <- c("Day89_Topic_14_top15", "Day94_Topic_27_top15","Day110_Topic14_top15", "Day115_Topic_12_top15", "Day125_Topic_33_top15")

for(i in 1:5){
  print(seuobj1_2_8 <- AddModuleScore(seuobj1_2_8, features = list("genes" = top15s[i]), assay = "RNA", name = top15names[i]))
}

for(j in top15namesnew){
  print(FeaturePlot(seuobj1_2_8, features = c(j), min.cutoff = "q1"))
}

DimPlot(seuobj1_2_8, reduction = "umap", group.by = "Main_cluster_name", label = TRUE)

DimPlot(seuobj1_2_8, reduction = "umap", group.by = "sub_cluster_name", label = TRUE)



seuobj1_2_8@meta.data$Day89_Topic_14_top151[is.na(seuobj1_2_8@meta.data$j)] <- 0
seuobj1_2_8@meta.data$Day94_Topic_27_top151[is.na(seuobj1_2_8@meta.data$j)] <- 0
seuobj1_2_8@meta.data$Day110_Topic14_top151[is.na(seuobj1_2_8@meta.data$j)] <- 0
seuobj1_2_8@meta.data$Day115_Topic_121[is.na(seuobj1_2_8@meta.data$j)] <- 0
seuobj1_2_8@meta.data$Day125_Topic_331[is.na(seuobj1_2_8@meta.data$j)] <- 0

DimPlot(seuobj1_2_8, group.by = "orig.ident")
colnames(obj@meta.data) <- gsub('old_column_name', 'new_column_name', colnames(obj@meta.data))

for(x in 1:5){
  print(colnames(seuobj1_2_8@meta.data) <- gsub(top15names[x], top15namesnew[x], colnames(seuobj1_2_8@meta.data)))
}



saveRDS(seuobj_all_days, file = "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/seuobj_all_days_kaustubh.RDS")



#pre-process integrated day 115 125 object












