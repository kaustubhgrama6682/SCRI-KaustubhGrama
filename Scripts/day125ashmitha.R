#load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)

#load data
seuobj125 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj125.RDS")





#astrocyte
seuobj125_astrocyte <- subset(seuobj125, subset = Main_cluster_name == 'Astrocytes')


seuobj125_astrocyte <- NormalizeData(seuobj125_astrocyte)
seuobj125_astrocyte <- FindVariableFeatures(seuobj125_astrocyte, selection.method = "vst", nfeatures = 2000)



#Scaling---------------------------------------------------
all.genes <- rownames(seuobj125_astrocyte)
seuobj125_astrocyte <- ScaleData(seuobj125_astrocyte, features = all.genes)

#PCA--------------------------------------------------------
seuobj125_astrocyte <- RunPCA(seuobj125_astrocyte, features = VariableFeatures(object = seuobj125_astrocyte))
DimHeatmap(seuobj125_astrocyte, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj125_astrocyte)




#granule neuron
seuobj125_granule_neuron <- subset(seuobj125, subset = Main_cluster_name == 'Granule neurons')

seuobj125_granule_neuron <- NormalizeData(seuobj125_granule_neuron)
seuobj125_granule_neuron <- FindVariableFeatures(seuobj125_granule_neuron, selection.method = "vst", nfeatures = 2000)



#Scaling---------------------------------------------------
all.genes <- rownames(seuobj125_granule_neuron)
seuobj125_granule_neuron <- ScaleData(seuobj125_granule_neuron, features = all.genes)

#PCA--------------------------------------------------------
seuobj125_granule_neuron <- RunPCA(seuobj125_granule_neuron, features = VariableFeatures(object = seuobj125_granule_neuron))
DimHeatmap(seuobj125_granule_neuron, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj125_granule_neuron)




#astrocyte and granule neuron


seuobj125_astrocyte_granule <- subset(seuobj125, subset = Main_cluster_name == "Astrocytes" | Main_cluster_name == "Granule neurons")

seuobj125_astrocyte_granule <- NormalizeData(seuobj125_astrocyte_granule)
seuobj125_astrocyte_granule <- FindVariableFeatures(seuobj125_astrocyte_granule, selection.method = "vst", nfeatures = 2000)



#Scaling---------------------------------------------------
all.genes <- rownames(seuobj125_astrocyte_granule)
seuobj125_astrocyte_granule <- ScaleData(seuobj125_astrocyte_granule, features = all.genes)

#PCA--------------------------------------------------------
seuobj125_astrocyte_granule <- RunPCA(seuobj125_astrocyte_granule, features = VariableFeatures(object = seuobj125_astrocyte_granule))
DimHeatmap(seuobj125_astrocyte_granule, dims = 1:18, cells = 500, balanced = TRUE)
ElbowPlot(seuobj125_astrocyte_granule)





#ADD MODULE SCORES

PC_1_Positive <- c("GRIK2", "RELN", "ERBB4", "LINC01036", "GPC6", "CNTN5", "KCND2", "LINGO2", "SYT1", "KCNIP4","CCSER1", "DLGAP1", "ROBO1", "CADPS2", "KLHL29", "CNTNAP2", "ZNF804B", "UNC5C", "ZNF385B", "XKR4","CNTN1", "MGAT4C", "CDH12", "MYRIP", "RP11-202H2.1", "THSD7B", "KCNH7", "DSCAML1", "LPHN2", "IL1RAPL1")
PC_1_Negative <-c("PTN", "SLC1A3", "NCKAP5", "SOX2-OT", "NPAS3", "TNC", "GRID2", "CTNNA3", "QKI", "PTPRM", "SOX6", "SLC35F1", "GPM6B", "PTPRZ1", "FAT3", "DMD", "EGFEM1P", "CDH4", "RP11-820L6.1", "LINC00478", "DCC", "MEGF10", "SMOC1", "SLC4A4", "PAX3", "ZBTB20", "PRKG1", "PTPRT", "GABRG3", "PLCE1")

PC_2_Positive <- c("LPPR1", "FSTL4", "IGFBP5", "NTN1", "INSM1", "TRHDE", "RNF220", "PDE1C", "ASPM", "HNRNPA1", "EGR1", "SFRP1", "EBF1", "PTPRK", "PCBP3", "PLXNB2", "RPS29", "EPHA7", "HS3ST3A1", "CYYR1", "CDK6", "CENPF", "EPHA5", "ROR2", "ZEB1", "UNC5C", "TOP2A", "RPS2", "RP11-154H12.3", "HS3ST5")
PC_2_Negative <- c("PCDH9", "CSMD3", "IL1RAPL1", "CTC-340A15.2", "ZNF804B", "TENM2", "LRRC4C", "CDH12", "LRRTM4", "EPHA6", "RP11-202H2.1", "MGAT4C", "KIRREL3", "KCNIP4", "NKAIN2", "HS3ST4", "KAZN", "DPP6", "DOK6", "SGCZ", "KCNH7", "LUZP2", "MDGA2", "CSMD1", "RP11-384F7.2", "CNTNAP5", "RGS6", "CDH18", "CTC-535M15.2", "CDH20")

seuobj125_astrocyte_granule <- AddModuleScore(seuobj125_astrocyte_granule, features = list("PC"=PC_1_Positive), assay = "RNA", name = "PC_1_Positive")
seuobj125_astrocyte_granule <- AddModuleScore(seuobj125_astrocyte_granule, features = list("PC"=PC_1_Negative), assay = "RNA", name = "PC_1_Negative")

seuobj125_astrocyte_granule <- AddModuleScore(seuobj125_astrocyte_granule, features = list("PC"=PC_2_Positive), assay = "RNA", name = "PC_2_Positive")
seuobj125_astrocyte_granule <- AddModuleScore(seuobj125_astrocyte_granule, features = list("PC"=PC_2_Negative), assay = "RNA", name = "PC_2_Negative")

pcs <- c("PC_1_Positive1", "PC_1_Negative1", "PC_2_Positive1", "PC_2_Negative1")

#Feature Plot
for(pc in pcs){
  print(FeaturePlot(seuobj125_astrocyte_granule, features = c(pc), min.cutoff = "q1"))
}

DimPlot(seuobj125_astrocyte_granule, group.by = "Main_cluster_name", label = TRUE)





















