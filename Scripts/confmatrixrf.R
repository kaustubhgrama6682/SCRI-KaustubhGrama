
seuobj94 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj94.RDS")
seuobj110 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj110.RDS")

astrocyte_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/astrocyte_rf.RDS")
II_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Inhibitory_interneuronrf.RDS")
UBC_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Unipolar_brush_cell_rf.RDS")
oligodendrocyte_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Oligodendrocyte_rf.RDS")
gn_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/gnrf.RDS")
vascular_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Vascular_endothelial_cell_rf.RDS")
microglia_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Microglia_rf.RDS")
purkinje_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Purkinje_neuron_rf.RDS")
SLC_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/SLC24A4_PEX5L_positive_cell_rf.RDS")

model_list <- list("Astrocytes" = astrocyte_model, 
                   "Inhibitory interneurons" = II_model, 
                   "Unipolar brush cells" = UBC_model, 
                   "Oligodendrocytes" = oligodendrocyte_model,
                   "Granule neurons" = gn_model,
                   "Vascular endothelial cells" = vascular_model,
                   "Microglia" = microglia_model,
                   "Purkinje neurons" = purkinje_model,
                   "SLC24A4_PEX5L positive cells" = SLC_model)



day94_gene_expression_data <- as.data.frame(seuobj94@assays$RNA$data)
day94_gene_expression_data <- mutate_all(day94_gene_expression_data, function(x) as.numeric(as.character(x)))
# switching cols and rows by transposing
day94_gene_expression_data <- t(day94_gene_expression_data)


# make data remain as data.frame
day94_gene_expression_data <- data.frame(day94_gene_expression_data)


# round vals to 2 spots
day94_gene_expression_data <- (round(day94_gene_expression_data, 2))


# add cell id col to metadata of seuobj
seuobj94$cellid <- rownames(seuobj94@meta.data)


# making a metadata dataframe of seuobj with main cluster and cell id
metadata_df <- seuobj94@meta.data[, c("Main_cluster_name", "cellid")]


# adding cell ids to metadata
day94_gene_expression_data$cellid <- rownames(day94_gene_expression_data)


# combining data and metadata
day94_gene_expression_data <- cbind(day94_gene_expression_data, metadata_df)


# removing the cell id from the expression dataframe
day94_gene_expression_data$cellid <- NULL





day110_gene_expression_data <- as.data.frame(seuobj110@assays$RNA$data)
day110_gene_expression_data <- mutate_all(day110_gene_expression_data, function(x) as.numeric(as.character(x)))
# switching cols and rows by transposing
day110_gene_expression_data <- t(day110_gene_expression_data)


# make data remain as data.frame
day110_gene_expression_data <- data.frame(day110_gene_expression_data)


# round vals to 2 spots
day110_gene_expression_data <- (round(day110_gene_expression_data, 2))


# add cell id col to metadata of seuobj
seuobj110$cellid <- rownames(seuobj110@meta.data)


# making a metadata dataframe of seuobj with main cluster and cell id
metadata_df <- seuobj110@meta.data[, c("Main_cluster_name", "cellid")]


# adding cell ids to metadata
day110_gene_expression_data$cellid <- rownames(day110_gene_expression_data)


# combining data and metadata
day110_gene_expression_data <- cbind(day110_gene_expression_data, metadata_df)


# removing the cell id from the expression dataframe
day110_gene_expression_data$cellid <- NULL



seuobj_115_125 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/seuobj_115_125.RDS")

#creating data frame of gene expression data
day115_125_gene_expression_data <- as.data.frame(seuobj_115_125@assays$integrated$data)

# switching cols and rows by transposing
day115_125_gene_expression_data <- t(day115_125_gene_expression_data)


# make data remain as data.frame
day115_125_gene_expression_data <- data.frame(day115_125_gene_expression_data)


# round vals to 2 spots
day115_125_gene_expression_data <- (round(day115_125_gene_expression_data, 2))


# add cell id col to metadata of seuobj
seuobj_115_125$cellid <- rownames(seuobj_115_125@meta.data)


# making a metadata dataframe of seuobj with main cluster and cell id
metadata_df <- seuobj_115_125@meta.data[, c("Main_cluster_name", "cellid")]


# adding cell ids to metadata
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)


# combining data and metadata
day115_125_gene_expression_data <- cbind(day115_125_gene_expression_data, metadata_df)


# removing the cell id from the expression dataframe
day115_125_gene_expression_data$cellid <- NULL

#making the training dat groping by main cluste name and filtering out
training.data <- day115_125_gene_expression_data %>%
  group_by(Main_cluster_name) %>%
  filter(row_number() <= 0.75 * n())






for(i in names(model_list)){
  
  day110_gene_expression_data$celltype <- 1
  day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == i] <- 0
  day94_gene_expression_data$celltype <- 1
  day94_gene_expression_data$celltype[day94_gene_expression_data$Main_cluster_name == i] <- 0
  
  probabilities110 <- predict(model_list[[i]], day110_gene_expression_data) 
  probabilities94 <- predict(model_list[[i]], day94_gene_expression_data) 
  
  original_odds <- 72266/217612
  undersample_odds <- sum(training.data$Main_cluster_name == i)/sum(training.data$Main_cluster_name != i)
  
  scoring_odds110 <- probabilities110/(1-probabilities110)
  scoring_odds94 <- probabilities94/(1-probabilities94)
  
  adjusted_odds110 <- scoring_odds110 * (original_odds/undersample_odds)
  adjusted_odds94 <- scoring_odds94 * (original_odds/undersample_odds)
  
  adjusted_probability110 = 1/(1+(1/adjusted_odds110))
  adjusted_probability94 = 1/(1+(1/adjusted_odds94))
  
  predicted.classes110 <- ifelse(adjusted_probability110 > 0.6, 1, 0) #adjusted prob to change
  predicted.classes94 <- ifelse(adjusted_probability94 > 0.6, 1, 0) #adjusted prob to change
  
  assign(paste("conf_matrix", i, "rf110", sep = ""), confusionMatrix(table(predicted.classes110, day110_gene_expression_data$celltype)))
  assign(paste("conf_matrix", i, "rf94", sep = ""), confusionMatrix(table(predicted.classes94, day94_gene_expression_data$celltype)))
  
  
  seuobj110$probabilities_neg110 <- 1-adjusted_probability110
  seuobj94$probabilities_neg94 <- 1-adjusted_probability94
  
  pt <- FeaturePlot(seuobj110, features = "probabilities_neg110") + ggtitle(paste(i , "rf110nocutoff"))
  png(paste0("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/FeaturePlots/", i, "rf110nocutoff.png"))
  print(pt)
  dev.off()
  
  
  pt1 <- FeaturePlot(seuobj110, features = "probabilities_neg110", min.cutoff = "q1") + ggtitle(paste(i , "rf110cutoff"))
  png(paste0("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/FeaturePlots/", i, "rf110cutoff.png"))
  print(pt1)
  dev.off()
  
  pt2 <- FeaturePlot(seuobj94, features = "probabilities_neg94") + ggtitle(paste(i , "rf94nocutoff"))
  png(paste0("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/FeaturePlots/", i, "rf94nocutoff.png"))
  print(pt2)
  dev.off()
  
  pt3 <- FeaturePlot(seuobj94, features = "probabilities_neg94", min.cutoff = "q1") + ggtitle(paste(i , "rf94cutoff"))
  png(paste0("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/FeaturePlots/", i, "rf94cutoff.png"))
  print(pt3)
  dev.off()
  
}















