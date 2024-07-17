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

# log plots
count <- unique(day115_125_gene_expression_data[,1])
count_df <- data.frame(count=count, percent=rep(NA, length(count)), logodds=rep(NA,length(count)))

# n_astro_list <- c()


rsq_df <- data.frame(gene = colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"], rsquare = rep(NA, 2000))

rsq_df[]


basepath <- "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/LogOdds_plots_Astrocytes/highest_r^2/"
for(gene in colnames(day115_125_gene_expression_data)[colnames(day115_125_gene_expression_data) != "Main_cluster_name"][1:100]){
  
  message(gene)

  count <- unique(day115_125_gene_expression_data[,gene])
  
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  message('Count df created')
  
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
    
    
    model <- lm(logodds ~ count, data = count_df)
    model_summary <- glance(model)
    r_squared <- model_summary$r.squared
    
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
  

rsq_df <- rsq_df[order(-rsq_df$rsquare),]

for(gene in rsq_df$gene[1:30]){
  message(gene)
  
  count <- unique(day115_125_gene_expression_data[,gene])
  
  count_df <- data.frame("count" = count, percent = rep(0, length(count)), logodds = rep(0, length(count)))
  
  message('Count df created')
  
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



    png(paste0(basepath, gene, ".png"), width = 300, height = 300)
    print(plot)
    dev.off()
  }
}


saveRDS(plot, file = "/Users/kaustubhgrama/Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/plot.RDS")
  




