
library(ROCR)
library(car)
seuobj110 <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/drive-download-20240702T163704Z-001/seuobj110.RDS")
astrocyte_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/model_lm_astrocyte.RDS")
II_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/model_lm_inhibitory_interneurons.RDS")
UBC_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/step.model.both_unipolar_brush_cells.RDS")
oligodendrocyte_model <- readRDS("/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/models/step.model.both_oligodendrocytes.RDS")


models <- list("Astrocytes" = astrocyte_model, 
               "Inhibitory interneurons" = II_model, 
               "Unipolar brush cells" = UBC_model, 
               "Oligodendrocytes" = oligodendrocyte_model)

models_str <- c("astrocyte", "II", "UBC", "oligodendrocyte")
cell_types <- c("Astrocytes", "Inhibitory interneurons", "Unipolar brush cells", "Oligodendrocytes")



#day 110 dat frame

#testing with day110 astrocyte model


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





#adjusted probability calc 




# 
# probabilities <- model_lm %>% predict(day110_gene_expression_data, type = "response") 
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



for(i in names(models)){
  
  
  
  probabilities <- models[[i]] %>% predict(day110_gene_expression_data, type = "response")
  undersample_odds <- sum(training.data$Main_cluster_name == i)/sum(training.data$Main_cluster_name != i)
  scoring_odds <- probabilities/(1-probabilities)
  adjusted_odds <- scoring_odds * (original_odds/undersample_odds)
  adjusted_probability = 1/(1+(1/adjusted_odds))
  colname <- paste(i, "probability")
  seuobj110@meta.data[, colname] <- adjusted_probability
  
  
  
}


df <- seuobj110@meta.data[,grep(" probability", colnames(seuobj110@meta.data))]
df <- colnames(df)[apply(df,1,which.max)]

seuobj110@meta.data$max_probability <- df 



