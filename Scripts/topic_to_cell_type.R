# Simple Gene Topic Analysis - Percentage of Marker Genes in Each Topic
# This calculates what percentage of cell type marker genes are found in each topic

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
# ==================== DATA LOADING ====================

# Load your topic genes from Excel
load_topic_genes <- function(file_path) {
  # Read the Excel file
  topic_data <- read_excel(file_path)
  
  # Convert to long format: each row = one gene-topic pair
  topic_genes <- topic_data %>%
    pivot_longer(cols = starts_with("Topic_"), 
                 names_to = "Topic", 
                 values_to = "Gene") %>%
    filter(!is.na(Gene)) %>%
    select(Topic, Gene)
  
  return(topic_genes)
}

# Your cerebellar markers
load_cerebellar_markers <- function() {
  celltypemarkers <- list(
    # Cell type specific lists
    rhombic_lip_progenitors = c("ATOH1", "PAX6", "LMX1A", "LHX9", "LHX2", "WLS",
                                "BMP7", "GDF7", "MSX1", "MSX2", "OLIG3"),
    granule_cell_precursors_proliferating = c("ATOH1", "PAX6", "MEIS1", "MEIS2",
                                              "CCND2", "MYCN", "SHH", "GLI1", "GLI2",
                                              "SMO", "PCNA", "KI67", "CDK6"),
    granule_cell_precursors_transitioning = c("NEUROD1", "TBR2", "NSCL1", "NSCL2",
                                              "BTG2", "GADD45A", "CDKN1B", "CDKN1C"),
    migrating_granule_cells = c("TAG1", "DCX", "ASTN1", "ASTN2", "SEMA6A", "PLEXINA2",
                                "UNC5C", "DCC", "EPHB2", "EFNB1", "EFNB2", "BDNF"),
    early_granule_neurons = c("NEUROD2", "CBLN1", "CBLN3", "EN1", "EN2", "ETV1",
                              "NFIB", "ZFP423", "TIAM1", "ARHGEF7"),
    mature_granule_neurons = c("GABRA6", "GABRD", "PDE1C", "CALB1", "SLC1A6",
                               "GRIN2C", "GRIA4", "KCNC3", "CACNG2", "HOMER3",
                               "SHANK2", "PCP2", "ITPR1"),
    ventricular_zone_progenitors = c("SOX2", "PAX3", "NESTIN", "HES1", "HES5",
                                     "NOTCH1", "DLL1", "JAG1", "RBPJ"),
    radial_glia = c("FABP7", "SLC1A3", "TNC", "VIM", "SOX2", "SOX9", "HOPX",
                    "PTPRZ1", "FGFR3", "EGFR"),
    astrocyte_precursors = c("SOX9", "NFIA", "NFIB", "ID3", "ID4", "FGFR3",
                             "OLIG2", "ASCL1", "COUP-TFI", "COUP-TFII"),
    bergmann_glia = c("GFAP", "S100B", "SLC1A3", "AQP4", "AMOT", "SEPT4",
                      "GLT1", "DIO2", "LRIG1", "FABP7", "SOX2"),
    velate_astrocytes = c("GFAP", "S100B", "KIR4.1", "KCNJ10", "GJA1", "GJB6",
                          "SLC1A3", "NTSR2", "PLA2G7"),
    fibrous_astrocytes = c("GFAP", "ALDH1L1", "AQP4", "S100B", "CD44", "LCN2",
                           "TIMP1", "SERPINA3N", "CXCL10"),
    # Main lineage lists
    rhombic_lip = c("ATOH1", "PAX6", "TBR2", "TBR1", "BARHL1", "ZIC1", "ZIC2",
                    "LMX1A", "LHX9", "LHX2", "MEIS1", "MEIS2", "WLS", "PTF1A"),
    granule_neuron = c("ATOH1", "PAX6", "BARHL1", "ZIC1", "ZIC2", "NEUROD1",
                       "GABRA6", "GABRD", "CBLN1", "CBLN3", "NSCL1", "NSCL2",
                       "EN1", "EN2", "ETV1", "UNC5C", "SEMA6A", "PLEXINA2",
                       "TAG1", "PDE1C", "CALB1", "SLC1A6", "GRIN2C"),
    astrocyte = c("SOX9", "NFIA", "NFIB", "ID3", "ID4", "HES1", "HES5",
                  "GFAP", "S100B", "SLC1A3", "AQP4", "ALDH1L1", "FABP7",
                  "VIM", "TNC", "SPARC", "CLU", "APOE", "GJA1", "GJB6"),
    # Biological processes
    proliferation = c("PCNA", "KI67", "CCND1", "CCND2", "CCNE1", "CDK2", "CDK4",
                      "CDK6", "MYCN", "MYC", "E2F1", "E2F3", "MCM2", "MCM6"),
    cell_cycle_exit = c("CDKN1A", "CDKN1B", "CDKN1C", "GADD45A", "BTG2", "TIS21",
                        "RB1", "RBL1", "RBL2", "GAS1"),
    differentiation = c("NEUROD1", "NEUROD2", "NFIA", "NFIB", "SOX9", "ASCL1",
                        "OLIG2", "MASH1", "NGN1", "NGN2", "POU3F2", "POU3F3"),
    migration = c("DCX", "ASTN1", "ASTN2", "TAG1", "PAFAH1B1", "RELN", "DAB1",
                  "CDK5", "CDK5R1", "NDEL1", "NDE1", "MARK1", "FYN"),
    axon_guidance = c("ROBO1", "ROBO2", "SLIT1", "SLIT2", "DCC", "UNC5C", "NETRIN1",
                      "SEMA3A", "SEMA3F", "SEMA6A", "PLEXINA1", "PLEXINA2", "NRP1",
                      "NRP2", "EPHB2", "EFNB1", "EFNB2"),
    synaptogenesis = c("CBLN1", "CBLN3", "NRXN1", "NRXN2", "NRXN3", "NLGN1",
                       "NLGN2", "NLGN3", "LRRTM1", "LRRTM2", "SYN1", "SYN2",
                       "SYP", "SNAP25", "VAMP2", "STX1A"),
    synaptic_maturation = c("ARC", "HOMER1", "HOMER3", "SHANK1", "SHANK2", "DLG4",
                            "CAMK2A", "CAMK2B", "GRIN2A", "GRIN2C", "GRIA1",
                            "GRIA4", "GRID2", "GRM1"),
    apoptosis_survival = c("BCL2", "BCL2L1", "BAX", "BAD", "CASP3", "CASP9",
                           "AKT1", "PTEN", "TRP53", "MDM2", "IGF1", "IGF1R",
                           "BDNF", "NTRK2"),
    metabolism = c("LDHA", "LDHB", "PFKFB3", "HK1", "HK2", "GLUT1", "GLUT3",
                   "MCT1", "MCT2", "COX4I1", "COX4I2", "NDUFA4", "ATP5A1"),
    epigenetic_regulation = c("EZH2", "BMI1", "RING1B", "SUZ12", "KDM6A", "KDM6B",
                              "HDAC1", "HDAC2", "HDAC3", "MLL1", "DOT1L", "DNMT1",
                              "DNMT3A", "DNMT3B"),
    transcriptional_regulation = c("REST", "COREST", "CTCF", "YY1", "SP1", "SP3",
                                   "CREB1", "ATF3", "JUN", "FOS", "MEF2A", "MEF2C",
                                   "MEF2D", "NFAT1", "NFATC4"),
    rna_processing = c("RBFOX1", "RBFOX2", "RBFOX3", "NOVA1", "NOVA2", "PTBP1",
                       "PTBP2", "SRRM4", "ELAVL2", "ELAVL3", "ELAVL4", "QKI"),
    calcium_signaling = c("CACNA1A", "CACNA1G", "CACNG2", "CALB1", "CALB2",
                          "PVALB", "ITPR1", "ITPR3", "RYR1", "RYR3", "ATP2B1",
                          "ATP2B2", "NCX1", "NCX2"),
    growth_factor_signaling = c("SHH", "SMO", "GLI1", "GLI2", "GLI3", "BDNF",
                                "NTF3", "NTRK2", "NTRK3", "IGF1", "IGF1R", "FGF8",
                                "FGFR1", "FGFR2", "FGFR3"),
    cell_adhesion = c("NCAM1", "NCAM2", "L1CAM", "CHL1", "NRCAM", "CADM1",
                      "CADM2", "CDH2", "CDH4", "CDH8", "CDH10", "PCDH10",
                      "CNTN1", "CNTN2"),
    extracellular_matrix = c("TNC", "TNR", "LAMA1", "LAMB1", "LAMC1", "COL4A1",
                             "NCAN", "BCAN", "VCAN", "ACAN", "HSPG2", "DAG1",
                             "ITGB1", "ITGA6"),
    glial_differentiation = c("OLIG1", "OLIG2", "SOX8", "SOX9", "SOX10", "NKX2.2",
                              "ID2", "ID3", "ID4", "STAT3", "SMAD1", "SMAD5",
                              "TCF7L2", "HEY1", "HEY2"),
    myelination = c("MBP", "PLP1", "MAG", "MOG", "CNP", "OPALIN", "CLAUDIN11",
                    "FA2H", "UGT8A", "GAL3ST1", "ASPA"),
    # Temporal expression patterns
    early_development = c("ATOH1", "PAX6", "LMX1A", "SOX2", "HES1", "HES5"),
    mid_development = c("NEUROD1", "DCX", "TAG1", "NFIA", "SOX9"),
    late_development = c("GABRA6", "CBLN1", "EN1", "EN2", "GFAP", "S100B"),
    maturation = c("PDE1C", "GRM1", "GRID2", "ALDH1L1", "AQP4")
  )
  
  return(celltypemarkers)
}

# ==================== MAIN ANALYSIS FUNCTION ====================

# Calculate percentage overlap between topics and cell type markers
calculate_percentage_overlap <- function(topic_genes_df, marker_list) {
  
  # Create results data frame
  results <- data.frame()
  
  # Get unique topics
  topics <- unique(topic_genes_df$Topic)
  
  # Calculate percentage overlap for each topic-celltype pair
  for (topic in topics) {
    # Get genes in this topic
    topic_gene_list <- topic_genes_df$Gene[topic_genes_df$Topic == topic]
    
    for (cell_type in names(marker_list)) {
      # Get marker genes for this cell type
      marker_gene_list <- marker_list[[cell_type]]
      
      # Calculate overlap
      overlapping_genes <- intersect(topic_gene_list, marker_gene_list)
      overlap_count <- length(overlapping_genes)
      total_markers <- length(marker_gene_list)
      percentage_overlap <- (overlap_count / total_markers) * 100
      
      # Store results
      results <- rbind(results, data.frame(
        Topic = topic,
        Cell_Type = cell_type,
        Overlap_Count = overlap_count,
        Total_Markers = total_markers,
        Percentage_Overlap = percentage_overlap,
        Overlapping_Genes = paste(overlapping_genes, collapse = ", ")
      ))
    }
  }
  
  return(results)
}

# ==================== SIMPLE ASSIGNMENT FUNCTION ====================

# Find best cell type assignment for each topic
assign_topics_simple <- function(percentage_results, min_percentage = 10) {
  
  # Find best assignment for each topic
  best_assignments <- percentage_results %>%
    filter(Percentage_Overlap >= min_percentage) %>%  # Only consider meaningful overlaps
    group_by(Topic) %>%
    slice_max(Percentage_Overlap, n = 1) %>%
    ungroup() %>%
    select(Topic, Cell_Type, Percentage_Overlap, Overlap_Count, Total_Markers, Overlapping_Genes)
  
  return(best_assignments)
}

# ==================== SIMPLE VISUALIZATION ====================

# Create simple heatmap of percentage overlaps
plot_percentage_heatmap <- function(percentage_results) {
  
  # Pivot to wide format for heatmap
  heatmap_data <- percentage_results %>%
    select(Topic, Cell_Type, Percentage_Overlap) %>%
    pivot_wider(names_from = Cell_Type, values_from = Percentage_Overlap) %>%
    column_to_rownames("Topic") %>%
    as.matrix()
  
  # Create heatmap
  library(pheatmap)
  pheatmap(heatmap_data,
           main = "Percentage of Marker Genes Found in Each Topic",
           color = colorRampPalette(c("white", "red"))(100),
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           display_numbers = TRUE,
           number_format = "%.1f",
           fontsize_number = 6,
           fontsize = 8)
}

# Create bar plot of best assignments
plot_best_assignments <- function(best_assignments) {
  library(ggplot2)
  
  ggplot(best_assignments, aes(x = reorder(Topic, Percentage_Overlap), 
                               y = Percentage_Overlap, 
                               fill = Cell_Type)) +
    geom_col() +
    coord_flip() +
    labs(title = "Best Cell Type Assignment for Each Topic",
         x = "Topic", 
         y = "Percentage of Marker Genes Found (%)",
         fill = "Assigned Cell Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# ==================== MAIN FUNCTION ====================

# Simple analysis pipeline
run_simple_analysis <- function(excel_file_path) {
  
  cat("Loading data...\n")
  # Load topic genes
  topic_genes <- load_topic_genes(excel_file_path)
  
  # Load marker genes
  marker_list <- load_cerebellar_markers()
  
  cat("Calculating percentage overlaps...\n")
  # Calculate percentage overlaps
  percentage_results <- calculate_percentage_overlap(topic_genes, marker_list)
  
  cat("Finding best assignments...\n")
  # Find best assignments
  best_assignments <- assign_topics_simple(percentage_results, min_percentage = 5)
  
  cat("Creating visualizations...\n")
  # Create plots
  plot_percentage_heatmap(percentage_results)
  plot_best_assignments(best_assignments)
  
  cat("Analysis complete!\n")
  
  return(list(
    percentage_results = percentage_results,
    best_assignments = best_assignments,
    topic_genes = topic_genes
  ))
}

# ==================== USAGE EXAMPLES ====================

cat("=== SIMPLE GENE TOPIC ANALYSIS ===\n\n")
cat("USAGE:\n")
cat("1. Install packages: install.packages(c('readxl', 'dplyr', 'tidyr', 'ggplot2', 'pheatmap'))\n")
cat("2. Run analysis: results <- run_simple_analysis('your_file.xlsx')\n")
cat("3. View results: View(results$best_assignments)\n")
cat("4. View all percentages: View(results$percentage_results)\n\n")

cat("EXAMPLE OUTPUT:\n")
cat("Topic_1  ->  granule_neuron (45.2% of markers found)\n")
cat("Topic_2  ->  astrocyte (32.1% of markers found)\n")
cat("Topic_3  ->  proliferation (67.8% of markers found)\n\n")

# Function to export simple results
export_simple_results <- function(results, output_file = "topic_assignments.csv") {
  
  # Export best assignments
  write.csv(results$best_assignments, output_file, row.names = FALSE)
  
  # Export full percentage matrix
  full_file <- gsub(".csv", "_full_percentages.csv", output_file)
  write.csv(results$percentage_results, full_file, row.names = FALSE)
  
  cat("Results exported to:", output_file, "and", full_file, "\n")
}

# Function to show top matches for each topic
show_top_matches <- function(percentage_results, top_n = 3) {
  
  top_matches <- percentage_results %>%
    group_by(Topic) %>%
    arrange(desc(Percentage_Overlap)) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  return(top_matches)
}

cat("ADDITIONAL FUNCTIONS:\n")
cat("# Export results\n")
cat("export_simple_results(results, 'my_assignments.csv')\n\n")
cat("# Show top 3 matches for each topic\n")

--------------------------------------------
  
  
  

cat("top_matches <- show_top_matches(results$percentage_results, top_n = 3)\n")
cat("View(top_matches)\n")







# Run the analysis
results <- run_simple_analysis("/Users/kaustubhgrama/Downloads/SCRI/TopGenes_Cao.xlsx")

# View best assignments for each topic
View(results$best_assignments)

# View all percentages
View(results$percentage_results)


days <- c("89", "94", "110", "115", "125")


for(day in days){
  results <- run_simple_analysis(paste0("/Users/kaustubhgrama/Downloads/SCRI/TopGenes_Cao", day,".xlsx"))
  saveRDS(results, file = paste0("/Users/kaustubhgrama/Downloads/SCRI/TopicResultsDay", day, ".RDS"))
  
}






