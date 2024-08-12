
library(ROCR)
library(car)
seuobj94 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj94.RDS")



astrocyte_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/astrocyte_rf.RDS")
II_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Inhibitory_interneuronrf.RDS")
UBC_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Unipolar_brush_cell_rf.RDS")
oligodendrocyte_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Oligodendrocyte_rf.RDS")
gn_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/gnrf.RDS")
vascular_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Vascular_endothelial_cell_rf.RDS")
microglia_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Microglia_rf.RDS")
purkinje_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/Purkinje_neuron_rf.RDS")
SLC_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/rf/drive-download-20240812T172845Z-001/SLC24A4_PEX5L_positive_cell_rf.RDS")

models_list <- list("Astrocytes" = astrocyte_model, 
               "Inhibitory interneurons" = II_model, 
               "Unipolar brush cells" = UBC_model, 
               "Oligodendrocytes" = oligodendrocyte_model,
               "Granule neurons" = gn_model,
               "Vascular endothelial cells" = vascular_model,
               "Microglia" = microglia_model,
               "Purkinje neurons" = purkinje_model,
               "SLC24A4_PEX5L positive cells" = SLC_model)






#testing with day94 astrocyte model


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





#adjusted probability calc 




# 
# probabilities <- model_lm %>% predict(day94_gene_expression_data, type = "response") 
# 
original_odds <- 72266/217612
# 
# undersample_odds <- sum(training.data$Main_cluster_name == "Astrocytes")/sum(training.data$Main_cluster_name != "Astrocytes")
# 
# scoring_odds <- probabilities/(1-probabilities)
# 
# adjusted_odds <- scoring_odds * (original_odds/undersample_odds)
# 
# adjusted_probability = 1/(1+(1/adjusted_odds))

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







for(i in names(models_list)){
  
  probabilities <-predict(models_list[[i]], day94_gene_expression_data)
  undersample_odds <- sum(training.data$Main_cluster_name == i)/sum(training.data$Main_cluster_name != i)
  scoring_odds <- probabilities/(1-probabilities)
  adjusted_odds <- scoring_odds * (original_odds/undersample_odds)
  adjusted_probability = 1/(1+(1/adjusted_odds))
  colname <- paste(i, "probability")
  seuobj94@meta.data[, colname] <- 1 - adjusted_probability
  
}


df <- seuobj94@meta.data[,grep(" probability", colnames(seuobj94@meta.data))]
df <- colnames(df)[apply(df,1,which.max)]

seuobj94@meta.data$max_probability <- df 


df3 <- seuobj94@meta.data[,grep("Day89_Topic_", colnames(seuobj94@meta.data))]
df3 <- colnames(df3)[apply(df3,1,which.max)]

seuobj94@meta.data$max_topic_day89 <- df3


df4 <- seuobj94@meta.data[,grep("Day94_Topic_", colnames(seuobj94@meta.data))]
df4 <- colnames(df4)[apply(df4,1,which.max)]

seuobj94@meta.data$max_topic_day94 <- df4


df5 <- seuobj94@meta.data[,grep("Topic__", colnames(seuobj94@meta.data))]
df5 <- colnames(df5)[apply(df5,1,which.max)]

seuobj94@meta.data$max_topic_day94 <- df5


df6 <- seuobj94@meta.data[,grep("Day115_Topic_", colnames(seuobj94@meta.data))]
df6 <- colnames(df6)[apply(df6,1,which.max)]

seuobj94@meta.data$max_topic_day115 <- df6


df7 <- seuobj94@meta.data[,grep("Day125_Topic_", colnames(seuobj94@meta.data))]
df7 <- colnames(df7)[apply(df7,1,which.max)]

seuobj94@meta.data$max_topic_day125 <- df7


plotdf <- table(seuobj94$Main_cluster_name, seuobj94$max_probability)

plotdf <- as.data.frame(plotdf)

ggplot(plotdf, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs predicted cluster name in Day 94")



plotdf <- table(seuobj94$Main_cluster_name, seuobj94$max_topic_day89)

plotdf <- as.data.frame(plotdf)
plotdf$Var2 <- gsub("Day89_", "", plotdf$Var2)
plotdf$Var2 <- factor(plotdf$Var2, levels = 
                        mixedsort(unique(plotdf$Var2)))

ggplot(plotdf, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 89 topics in day 94 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))




plotdf2 <- table(seuobj94$Main_cluster_name, seuobj94$max_topic_day94)

plotdf2 <- as.data.frame(plotdf2)
plotdf2$Var2 <- gsub("Day94_", "", plotdf2$Var2)
plotdf2$Var2 <- factor(plotdf2$Var2, levels = 
                         mixedsort(unique(plotdf2$Var2)))

ggplot(plotdf2, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 94 topics in day 94 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))


plotdf3 <- table(seuobj94$Main_cluster_name, seuobj94$max_topic_day115)

plotdf3 <- as.data.frame(plotdf3)
plotdf3$Var2 <- gsub("Day115_", "", plotdf3$Var2)
plotdf3$Var2 <- factor(plotdf3$Var2, levels = 
                         mixedsort(unique(plotdf3$Var2)))

ggplot(plotdf3, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 115 topics in day 94 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))



plotdf4 <- table(seuobj94$Main_cluster_name, seuobj94$max_topic_day125)

plotdf4 <- as.data.frame(plotdf4)
plotdf4$Var2 <- gsub("Day125_", "", plotdf4$Var2)
plotdf4$Var2 <- factor(plotdf4$Var2, levels = 
                         mixedsort(unique(plotdf4$Var2)))

ggplot(plotdf4, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Main cluster name vs day 125 topics in day 94 data") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1))









