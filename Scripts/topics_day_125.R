library(Seurat)

seuobj125 <- readRDS(file = "/Users/kaustubhgrama/Downloads/SCRI/seuobj125.RDS")

day89df <- as_data_frame(seuobj125@meta.data)

topic_cols_89 <- grep("^Day89_Topic_", names(day89df), value = TRUE)

# Create the three new columns
day89df$highest_topic <- apply(day89df[topic_cols_89], 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})

day89df$second_highest_topic <- apply(day89df[topic_cols_89], 1, function(x) {
  names(sort(x, decreasing = TRUE))[2]
})

day89df$third_highest_topic <- apply(day89df[topic_cols_89], 1, function(x) {
  names(sort(x, decreasing = TRUE))[3]
})



#94 topic in 89 object
day94df <-  as_data_frame(seuobj125@meta.data)

topic_cols_94 <- grep("^Day94_Topic_", names(day94df), value = TRUE)

day94df$highest_topic <- apply(day94df[topic_cols_94], 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})

day94df$second_highest_topic <- apply(day94df[topic_cols_94], 1, function(x) {
  names(sort(x, decreasing = TRUE))[2]
})

day94df$third_highest_topic <- apply(day94df[topic_cols_94], 1, function(x) {
  names(sort(x, decreasing = TRUE))[3]
})



#topic 110 in 89 object
day110df <-  as_data_frame(seuobj125@meta.data)

topic_cols_110 <- grep("^Day110_Topic_", names(day110df), value = TRUE)

day110df$highest_topic <- apply(day110df[topic_cols_110], 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})

day110df$second_highest_topic <- apply(day110df[topic_cols_110], 1, function(x) {
  names(sort(x, decreasing = TRUE))[2]
})

day110df$third_highest_topic <- apply(day110df[topic_cols_110], 1, function(x) {
  names(sort(x, decreasing = TRUE))[3]
})


#topic 115 in 89 object
day115df <-  as_data_frame(seuobj125@meta.data)

topic_cols_115 <- grep("^Day115_Topic_", names(day115df), value = TRUE)

day115df$highest_topic <- apply(day115df[topic_cols_115], 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})

day115df$second_highest_topic <- apply(day115df[topic_cols_115], 1, function(x) {
  names(sort(x, decreasing = TRUE))[2]
})

day115df$third_highest_topic <- apply(day115df[topic_cols_115], 1, function(x) {
  names(sort(x, decreasing = TRUE))[3]
})


#topic 125 in 89 object
day125df <-  as_data_frame(seuobj125@meta.data)

topic_cols_125 <- grep("^Topic_", names(day115df), value = TRUE)

day125df$highest_topic <- apply(day125df[topic_cols_125], 1, function(x) {
  names(sort(x, decreasing = TRUE))[1]
})

day125df$second_highest_topic <- apply(day125df[topic_cols_125], 1, function(x) {
  names(sort(x, decreasing = TRUE))[2]
})

day125df$third_highest_topic <- apply(day125df[topic_cols_125], 1, function(x) {
  names(sort(x, decreasing = TRUE))[3]
})


# Add the new columns back to the Seurat object metadata
seuobj125@meta.data$highest_topic <- day125df$highest_topic
seuobj125@meta.data$second_highest_topic <- day125df$second_highest_topic
seuobj125@meta.data$third_highest_topic <- day125df$third_highest_topic

# Subset to "Granule Neurons"
astrocyte_obj <- subset(seuobj125, subset = Main_cluster_name == "Astrocytes")

# Create UMAP DimPlots for each new column
# Plot 1: Highest topic
p1 <- DimPlot(astrocyte_obj, 
              reduction = "umap", 
              group.by = "highest_topic",
              label = TRUE,
              label.size = 3) +
  ggtitle("Highest Topic per Cell") +
  NoLegend()

# Plot 2: Second highest topic  
p2 <- DimPlot(astrocyte_obj, 
              reduction = "umap", 
              group.by = "second_highest_topic",
              label = TRUE,
              label.size = 3) +
  ggtitle("Second Highest Topic per Cell") +
  NoLegend()

# Plot 3: Third highest topic
p3 <- DimPlot(astrocyte_obj, 
              reduction = "umap", 
              group.by = "third_highest_topic",
              label = TRUE,
              label.size = 3) +
  ggtitle("Third Highest Topic per Cell") +
  NoLegend()

# Display plots individually
print(p1)
print(p2)
print(p3)





