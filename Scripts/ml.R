#figure out what happpebned heree
# #creating data frame of gene expression data
# day115_125_gene_expression_data <- as.data.frame(seuobj_115_125@assays$integrated$data)
# #switching the columns and the rows
# day115_125_gene_expression_data <- t(day115_125_gene_expression_data)
# #making sure it stays as a data frame
# day115_125_gene_expression_data <- data.frame(day115_125_gene_expression_data)
# #rounding the values in the data frame
# day115_125_gene_expression_data <- (round(day115_125_gene_expression_data, 2))
# #adding the cell id column to the meta data of the original seuobj
# seuobj_115_125$cellid <- rownames(seuobj_115_125@meta.data)
# #making a metadata data frame from the meta data of the original seuobj with the main cluster and cell id coluns
# metadata_df <- seuobj_115_125@meta.data[, c("Main_cluster_name", "cellid")]
# #adding the actual cell ids into the cellid column
# day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)
# #joining the metadata of the original seuobj and the new dataframe
# day115_125_gene_expression_data <- cbind(day115_125_gene_expression_data, metadata_df)
# #removing the cell id from the expression dataframE
# day115_125_gene_expression_data$cellid <- NULL
# 
# #log lpots
# count <- unique(day115_125_gene_expression_data[,1])
# count_df <- data.frame(count=count, percent=rep(NA, length(count)), logodds=rep(NA,length(count)))
# 
# for(n in unique(count_df$count)){
#   # count for all astrocytes
#   n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",1]==n)
#   print(n_astrocytes)
# }
# 
# #n_astro_list <- c()
# for(n in unique(count_df$count)[1]){
#   # count for all astrocytes
#   n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",1]==n)
#   n_astro_list <- c(n_astro_list, n_astrocytes)
# }
# table(n_astro_list) 
# 
# 
# for(n in count_df$count){
# 	# count for all astrocytes
# 	n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",1]==n)
# 	# probability of a cell being an astrocyte
# 	# number of astroctyes w uniq count value for spp1 / all cells w uniq count value for spp1
# 	percent <- n_astrocytes/ sum(day115_125_gene_expression_data[,1]==n)
# 	# log odds
# 	logodds <- log(percent/1-percent)
# 	# adding all variables to count_df
# 	count_df[count_df$count == n, "percent"] <- percent
# 	count_df[count_df$count == n, "logodds"] <- logodds
# }


library(ggplot2)
library(ggpmisc)
library(nlme)
library(broom)
library(dplyr)
library(car)
library(Seurat)

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

# # log plots
# count <- unique(day115_125_gene_expression_data[,1])
# count_df <- data.frame(count=count, percent=rep(NA, length(count)), logodds=rep(NA,length(count)))

# n_astro_list <- c()

#create a data frame to store r^2 values for each gene
rsq_df <- data.frame(gene = colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"], rsquare = rep(NA, 2000))

rsq_df[]

#setting basepath to save png files
basepath <- "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/LogOdds_plots_Astrocytes/highest_r^2/above0.2/"

#for loop to populate rsq_df
for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"]){
  #tells name of current gene
  message(gene)
  
  #filling count with gene expression values for all cells for this gene
  count <- unique(day115_125_gene_expression_data[,gene])
  
  #creating a dataframe with percent and logodds columns corresponding to probability that cell is an astrocyte depending on counts for that gene
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  #displays that coutn df was created
  message('Count df created')
  
  #for loop iterating through each unique count value
  for(n in unique(count_df$count)){
    
    # count for all astrocytes
    n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",gene]==n)
    
    # probability of a cell being an astrocyte
    # number of astrocytes w unique count value for spp1 / all cells
    percent <- n_astrocytes/ sum(day115_125_gene_expression_data[,gene]==n)
    
    # log odds
    logodds <- log(percent/(1-percent))
    
    # adding all variables to count_df
    count_df[count_df$count == n, "percent"] <- percent
    count_df[count_df$count == n, "logodds"] <- logodds
  }

  
  #removing all -Inf logodds
  count_df <- count_df[count_df$logodds != -Inf, ]
  #removing all Inf logodds
  count_df <- count_df[count_df$logodds != Inf, ]
  
  #only executing if there are more than 20 unique counts for the gene
  if(nrow(count_df) >= 20){
    
    #creating a linear model to model relation between count and logodds
    model <- lm(logodds ~ count, data = count_df)
    #getting summary of the model
    model_summary <- glance(model)
    #extracting r squared from the summary
    r_squared <- model_summary$r.squared
    
    #adding r62 value from the model into rsq_df
    rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
    
    # plot <- ggplot(count_df, aes(x = count, y = logodds)) +
    #   geom_point(stat = "identity") + geom_smooth(method = "lm") +
    #   ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
    # 
    # 
    # 
    # png(paste0(basepath, gene, ".png"), width = 300, height = 300)
    # print(plot)
    # dev.off()
  }
}
  
#ordering rsq_df by descending r^2 
rsq_df <- rsq_df[order(-rsq_df$rsquare),]

#iterating through genes with highest 30 r^2
for(gene in rsq_df$gene[1:30]){
  #getting the index of the gene
  row_index <- which(rsq_df$gene == gene)
  #only move forward if the r^2 value for this gene is above 0.2
  if(rsq_df$rsquare[row_index] >= 0.2){
    #say the gene name
    message(gene)
    
    #getting all the unique count value for that gene
    count <- unique(day115_125_gene_expression_data[,gene])
    
    
    #making a count data frame with columns for percent and logodds
    count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
    
    #saying that count df was created
    message('Count df created')
    
    #iterating through all the unique count values
    for(n in unique(count_df$count)){
      # count for all astrocytes
      n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",gene]==n)
      
      # probability of a cell being an astrocyte
      # number of astrocytes w unique count value for spp1 / all cells
      percent <- n_astrocytes/ sum(day115_125_gene_expression_data[,gene]==n)
      
      # log odds
      logodds <- log(percent/(1-percent))
      
      # adding all variables to count_df
      count_df[count_df$count == n, "percent"] <- percent
      count_df[count_df$count == n, "logodds"] <- logodds
    }
    # print(count_df)
    # table(count_df$logodds==-Inf)
    
    #removing -inf logodds values
    count_df <- count_df[count_df$logodds != -Inf, ]
    #removing inf logodds values
    count_df <- count_df[count_df$logodds != Inf, ]
    
    #if the number of unique count values is at least 20 then proceed
    if(nrow(count_df) >= 20){
      
      
      # model <- lm(logodds ~ count, data = count_df)
      # model_summary <- glance(model)
      # r_squared <- model_summary$r.squared
      # 
      # rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared
      
      #create a ggplot for the logodds and counts
      plot <- ggplot(count_df, aes(x = count, y = logodds)) +
        geom_point(stat = "identity") + geom_smooth(method = "lm") +
        ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
      
      
      
      png(paste0(basepath, gene, ".png"), width = 300, height = 300)
      print(plot)
      dev.off()
  }
  
  }
}

#save rds file
saveRDS(plot, file = "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/plot.RDS")
  
#iterating through the 30 rows in rsq_df
for(i in 1 : 30) {

  #if the r^2 value is above 0.2 proceed
  if (rsq_df$rsquare[i] >= 0.4) {
    
    print(paste0("COUNT",i))
    
    #getting the gene associated with that row
    gene <- rsq_df$gene[i]
    message(gene)

    #making a counts that have unique count value
    count <- unique(day115_125_gene_expression_data[,gene])
    #count df with logodds and percent columns
    count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))

    message('Count df created')

    #iterating through all unique count values
    for(n in unique(count_df$count)){
      # count for all astrocytes
      n_astrocytes <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Astrocytes",gene]==n)

      # probability of a cell being an astrocyte
      # number of astrocytes w unique count value for spp1 / all cells
      percent <- n_astrocytes/ sum(day115_125_gene_expression_data[,gene]==n)

      # log odds
      logodds <- log(percent/(1-percent))

      # adding all variables to count_df
      count_df[count_df$count == n, "percent"] <- percent
      count_df[count_df$count == n, "logodds"] <- logodds
    }
    # print(count_df)
    # table(count_df$logodds==-Inf)

    count_df <- count_df[count_df$logodds != -Inf, ]
    count_df <- count_df[count_df$logodds != Inf, ]

    if(nrow(count_df) >= 20){


      # model <- lm(logodds ~ count, data = count_df)
      # model_summary <- glance(model)
      # r_squared <- model_summary$r.squared
      #
      # rsq_df[rsq_df$gene==gene, "rsquare"] <- r_squared

      plot <- ggplot(count_df, aes(x = count, y = logodds)) +
        geom_point(stat = "identity") + geom_smooth(method = "lm") +
        ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))



      png(paste0(basepath, gene, "TOP30",".png"), width = 300, height = 300)
      print(plot)
      dev.off()
    }
  }

}

#astrocyte model
#getting the top genes from rsq_)df
topgenes_logodds <- c()
for(i in 1 : 30) {
  if (rsq_df$rsquare[i] >= 0.4) {
    topgenes_logodds <- append(topgenes_logodds, rsq_df$gene[i])
  }
}

astrocyte.markers <- c("S100B","GFAP", "APOE", "AQP4")
topgenes_logodds <- append(topgenes_logodds, astrocyte.markers)



#load dplyr
library(dplyr)

#making the cell id columns the same as the rownames
day115_125_gene_expression_data$cellid <- rownames(day115_125_gene_expression_data)

#making the training dat groping by main cluste name and filtering out
training.data <- day115_125_gene_expression_data %>%
  group_by(Main_cluster_name) %>%
  filter(row_number() <= 0.75 * n())

#making testing data from everything that isnt training
testing.data <- day115_125_gene_expression_data[!day115_125_gene_expression_data$cellid %in% training.data$cellid, ]


rownames(training.data) <- training.data$cellid
rownames(testing.data) <- testing.data$cellid


training.data$cellid <- NULL
testing.data$cellid <- NULL


#creating top genes from each
training.data.topgenes <- training.data[, c(topgenes_logodds, "Main_cluster_name")]
testing.data.topgenes <- testing.data[, c(topgenes_logodds, "Main_cluster_name")]

#creating binary predcictor column (astrocyte = 1, other = 0)
training.data.topgenes$celltype <- 0
training.data.topgenes$celltype[training.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1
 
testing.data.topgenes$celltype <- 0
testing.data.topgenes$celltype[testing.data.topgenes$Main_cluster_name == "Astrocytes"] <- 1




training.data.topgenes$Main_cluster_name <- NULL
testing.data.topgenes$Main_cluster_name <- NULL


#training glm with 0.75 GE data
model_lm <- glm(celltype ~. , data = training.data.topgenes, family = "binomial")
  
#testing GLM with 0.25 GE data
predictions <- model_lm %>% predict(testing.data.topgenes)
  
#determining the model performance
performance <- data.frame(RMSE = RMSE(predictions, testing.data.topgenes$celltype),
                          R2 = R2(predictions, testing.data.topgenes$celltype))


#view model summmary and performance
summary(model_lm)
performance

#vif to determine multicollinearity
car::vif(model_lm)

#stepwise AIC 
#stepwise regression --forward
library(MASS)


step.model.forward <- stepAIC(model_lm, direction = "forward", trace = 5)
summary(step.model.forward)

step.model.backward <- stepAIC(model_lm, direction = "backward", trace = 5)
summary(step.model.backward)

step.model.both <- stepAIC(model_lm, direction = "both", trace = 5)
summary(step.model.both)




#workflow for neurons


rsq_df_neurons <- data.frame(gene = colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"], rsquare = rep(NA, 2000))

rsq_df_neurons[]





for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"]){
  #tells name of current gene
  message(gene)
  #filling count with gene expression values for all cells for this gene
  count <- unique(day115_125_gene_expression_data[,gene])
  #creating a dataframe with percent and logodds columns corresponding to probability that cell is a granule neuron or inhibitory interneuron depending on counts for that gene
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  #displays that coutn df was created
  message('Count df created')
  
  #for loop iterating through each unique count value
  for(n in unique(count_df$count)){
    # count for all neurons
    n_ii <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Inhibitory interneurons",gene]==n)
    n_gn <- sum(day115_125_gene_expression_data[day115_125_gene_expression_data$Main_cluster_name == "Granule neurons",gene]==n)
    
    n_neurons <- n_ii + n_gn
    
    # probability of a cell being an neuron
    # number of neurons w unique count value for spp1 / all cells
    percent <- n_neurons/ sum(day115_125_gene_expression_data[,gene]==n)
    
    # log odds
    logodds <- log(percent/(1-percent))
    
    # adding all variables to count_df
    count_df[count_df$count == n, "percent"] <- percent
    count_df[count_df$count == n, "logodds"] <- logodds
  }
  # print(count_df)
  # table(count_df$logodds==-Inf)
  
  #removing all -Inf logodds
  count_df <- count_df[count_df$logodds != -Inf, ]
  #removing all Inf logodds
  count_df <- count_df[count_df$logodds != Inf, ]
  
  #only executing if there are more than 20 unique counts for the gene
  if(nrow(count_df) >= 20){
    
    #creating a linear model to model relation between count and logodds
    model <- lm(logodds ~ count, data = count_df)
    #getting summary of the model
    model_summary <- glance(model)
    #extracting r squared from the summary
    r_squared <- model_summary$r.squared
    
    #adding r62 value from the model into rsq_df
    rsq_df_neurons[rsq_df_neurons$gene==gene, "rsquare"] <- r_squared
    
    # plot <- ggplot(count_df, aes(x = count, y = logodds)) +
    #   geom_point(stat = "identity") + geom_smooth(method = "lm") +
    #   ggtitle(paste0("Log odds for ", gene, "\n(Astrocytes)")) + stat_poly_eq(use_label(c("eq", "R2")))
    # 
    # 
    # 
    # png(paste0(basepath, gene, ".png"), width = 300, height = 300)
    # print(plot)
    # dev.off()
  }
}

rsq_df_neurons <- rsq_df_neurons[order(-rsq_df_neurons$rsquare),]

topgenes_logodds_neurons <- c()
for(i in 1 : 30) {
  if (rsq_df_neurons$rsquare[i] >= 0.4) {
    topgenes_logodds_neurons <- append(topgenes_logodds_neurons, rsq_df_neurons$gene[i])
  }
}


training.data.topgenes.neurons <- training.data[, c(topgenes_logodds_neurons, "Main_cluster_name")]
testing.data.topgenes.neurons <- testing.data[, c(topgenes_logodds_neurons, "Main_cluster_name")]


training.data.topgenes.neurons$celltype <- 0
training.data.topgenes.neurons$celltype[training.data.topgenes.neurons$Main_cluster_name == "Granule neurons"] <- 1
training.data.topgenes.neurons$celltype[training.data.topgenes.neurons$Main_cluster_name == "Inhibitory interneurons"] <- 1


testing.data.topgenes.neurons$celltype <- 0
testing.data.topgenes.neurons$celltype[testing.data.topgenes.neurons$Main_cluster_name == "Granule neurons"] <- 1
testing.data.topgenes.neurons$celltype[testing.data.topgenes.neurons$Main_cluster_name == "Inhibitory interneurons"] <- 1


training.data.topgenes.neurons$Main_cluster_name <- NULL
testing.data.topgenes.neurons$Main_cluster_name <- NULL




#training glm with 0.75 GE data
model_lm_neurons <- glm(celltype ~. , data = training.data.topgenes.neurons, family = "binomial")

#testing GLM with 0.25 GE data
predictions.neurons <- model_lm_neurons %>% predict(testing.data.topgenes.neurons)

#determining the model performance
performance.neurons <- data.frame(RMSE = RMSE(predictions.neurons, testing.data.topgenes.neurons$celltype),
                          R2 = R2(predictions.neurons, testing.data.topgenes.neurons$celltype))


#view model summmary and performance
summary(model_lm_neurons)
performance.neurons

#vif to determine multicollinearity
car::vif(model_lm_neurons)

#stepwise AIC 
#stepwise regression --forward
library(MASS)


step.model.forward.neurons <- stepAIC(model_lm_neurons, direction = "forward", trace = 5)
summary(step.model.forward.neurons)

step.model.backward.neurons <- stepAIC(model_lm_neurons, direction = "backward", trace = 5)
summary(step.model.backward.neurons)

step.model.both.neurons <- stepAIC(model_lm_neurons, direction = "both", trace = 5)
summary(step.model.both.neurons)


day110_gene_expression_data$celltype <- 0
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Granule neurons"] <- 1
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Inhibitory interneurons"] <- 1

#day 110 predictions neuron model
# Making predictions
probabilities.neurons <- model_lm_neurons %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes.neurons <- ifelse(probabilities > 0.6, 1, 0)


# determining accuracy
table(predicted.classes.neurons == day110_gene_expression_data$celltype)
confusionMatrix(table(predicted.classes.neurons, day110_gene_expression_data$celltype))














#testing with day110 astrocyte model


day110_gene_expression_data <- as.data.frame(seuobj110@assays$RNA$data)

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


day110_gene_expression_data$celltype <- 0
day110_gene_expression_data$celltype[day110_gene_expression_data$Main_cluster_name == "Astrocytes"] <- 1

# Making predictions
probabilities <- model_lm %>% predict(day110_gene_expression_data, type = "response") 
predicted.classes <- ifelse(probabilities > 0.6, 1, 0)


# determining accuracy
table(predicted.classes == day110_gene_expression_data$celltype)
confusionMatrix(table(predicted.classes, day110_gene_expression_data$celltype))












