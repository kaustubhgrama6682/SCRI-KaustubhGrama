library(Seurat)
library(ggplot2)


seuobj125_astrocyte <- readRDS(file = "/Users/kaustubhgrama/Downloads/SCRI/seuobj125_astrocyte.RDS")

days <- c("89", "94", "110", "115", "125")

for(day in days){
  
  df <- as.data.frame(seuobj125_astrocyte@meta.data)
  
  topiccolname <- paste0("^Day", day, "_Topic_")
  
  topic_cols <- grep(topiccolname, names(df), value = TRUE)
  
  # Create the three new columns
  code <- paste0(
  "df$highest_topic", day, " <- apply(df[topic_cols], 1, function(x) {
    names(sort(x, decreasing = TRUE))[1]
  })
  
  df$second_highest_topic", day, "  <- apply(df[topic_cols], 1, function(x) {
    names(sort(x, decreasing = TRUE))[2]
  })
  
  df$third_highest_topic", day, "  <- apply(df[topic_cols], 1, function(x) {
    names(sort(x, decreasing = TRUE))[3]
  })"
  )
  
  eval(parse(text = code))
  
  code2 <- paste0(
  "seuobj125_astrocyte@meta.data$highest_topic", day, "  <- df$highest_topic", day, " 
  seuobj125_astrocyte@meta.data$second_highest_topic", day, "  <- df$second_highest_topic", day, " 
  seuobj125_astrocyte@meta.data$third_highest_topic", day, "  <- df$third_highest_topic", day
  )
  
  eval(parse(text = code2))
  
  
  
  
  
  
  
}
















pdf("/Users/kaustubhgrama/Downloads/SCRI/seuobj125_astrocyte_topic_UMAPs/seuobj125_astrocyte_topic_UMAPs.pdf")

for(day in days){
  code3 <- paste0(
    "p1 <- DimPlot(seuobj125_astrocyte, 
                reduction = 'umap', 
                group.by = 'highest_topic", day, "',
                label = TRUE,
                label.size = 3) +
    ggtitle('Highest Topic per Cell ", day, "') +
    NoLegend()

  p2 <- DimPlot(seuobj125_astrocyte, 
                reduction = 'umap', 
                group.by = 'second_highest_topic", day, "',
                label = TRUE,
                label.size = 3) +
    ggtitle('Second Highest Topic per Cell ", day, "') +
    NoLegend()
  

  p3 <- DimPlot(seuobj125_astrocyte, 
                reduction = 'umap', 
                group.by = 'third_highest_topic", day, "',
                label = TRUE,
                label.size = 3) +
    ggtitle('Third Highest Topic per Cell ", day, "') +
    NoLegend()"
  )
  
  eval(parse(text = code3))
  
 
  
  
  # basepath <- '/Users/kaustubhgrama/Downloads/SCRI/seuobj125_astrocyte_topic_UMAPs/'
  # 
  # saveRDS(p1, file = paste0(basepath, 'highest_topic', day, ".png"))
  # saveRDS(p2, file = paste0(basepath, 'second_highest_topic', day, ".png"))
  # saveRDS(p3, file = paste0(basepath, 'third_highest_topic', day, ".png"))
  
  print(p1)
  print(p2)
  print(p3)
  
}

dev.off()









