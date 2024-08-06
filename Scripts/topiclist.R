# MB_Topics <- MB_Topics[-1, ]
# 
# 
# 
# 
# for(i in 1:ncol(MB_Topics)){
#   variableName <- paste("topic",i, sep = "")
#   topics <- append(topics, assign(variableName, MB_Topics[, i]))
#   
# }
# 
# topics_df <- data.frame(topic1, topic2, topic3, topic4, topic5, topic6, topic7, topic8, topic9,
#             topic10, topic11, topic12, topic13, topic14, topic15, topic16, topic17, topic18, topic19,
#             topic20, topic21, topic22, topic23, topic24, topic25, topic26, topic27, topic28, topic29, 
#             topic30, topic31, topic32, topic33, topic34, topic35, topic36, topic37, topic38, topic39,
#             topic40, topic41, topic42, topic43, topic44, topic45, topic46, topic47, topic48, topic49, topic50)
# 
# 
# wnt_mutation <- c("CTNNB1" , "DDX3X", "SMARCA4", "CREBBP", "KMT2D")
# wntA <- c("TOP2A", "CDK1", "RRM2")
# wntC <- c("STMN2","KIF5C", "SYT11")
# wntB <- c("NME2", "HK2", "PGM5")
# wntD <- c("LRP4", "APCDD1")
# wntInG3 <- c("NFATC4", "GALNT14")
# paper_genes <- c(wnt_mutation, wntA, wntC, wntB, wntD, wntInG3)
# 

list_of_topics <- list()

for(i in 1:min(50, ncol(df))){
  column_list <- as.list(MB_Topics[, i])
  list_of_topics[[paste0("topic", i)]] <- column_list
}

# Print the names of the lists created
print(names(list_of_topics))

# each topic has list of genes

paper_genes <- c(
  "CTNNB1", "DDX3X", "SMARCA4", "CREBBP", "KMT2D",
  "MYC", "GABRA5", "GFI1", "GFI1B", "OTX2",
  "BMI1", "SOX2", "NFATC4", "GALNT14", "APC",
  "TP53", "TOP2A", "CDK1", "RRM2", "STMN2",
  "KIF5C", "SYT11", "NME2", "HK2", "PGM5",
  "LRP4", "APCDD1", "JUNB", "EGR1", "DKK2",
  "AXIN2", "WIF1", "ZIC1", "PAX6", "OLIG3",
  "NKD1", "BARHL1", "PDE1C", "PCSK9", "PTF1A",
  "NEUROG1", "ASCL1", "FOXD3", "BRN3A", "CBFA2T2",
  "CBFA2T3", "PRDM6", "UTX"
)



gene_info_df <- read_excel("Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/gene_info_df.xlsx")

gene_info_df$topics <- NA



# Iterate through each mutation
for (i in 1:length(paper_genes)) {
  
  mutation <- paper_genes[i]
  
  # Check in which topics the current mutation is present
  is_in_list <- sapply(list_of_topics, function(topic) mutation %in% topic)
  
  # Find the indices of topics containing the mutation
  topics_with_mutation <- which(is_in_list)
  
  # Print the result for the current mutation
  if (length(topics_with_mutation) > 0) {
    topic_names <- paste("Topic", topics_with_mutation, collapse = ", ")
    print(paste("Gene", mutation, "is in", topic_names))
    
    
    topics_with_mutation_string <- paste(topics_with_mutation, collapse = " ")
    
    gene_info_df$topics[i] <- topics_with_mutation_string
    
    
    
  } else {
    print(paste("Gene", mutation, "is not in any topic."))
    
    
    
    
  }
  

}

write_xlsx(gene_info_df, "Desktop/Computer_Science/R/Data/fetal_cerebellar_scData/gene_info_df.xlsx")






