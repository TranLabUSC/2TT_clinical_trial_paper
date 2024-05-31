library(ggplot2)
library(survival)
library(Seurat)
library(ggfortify)
library("Rcpp")
library(dplyr)
library(tidyr)

create_survival_data <- function(gmt_file, pathway_name, seurat_obj, survival_data) {
  # Read GMT file and extract genes for the specified pathway
  gmt <- readLines(gmt_file)
  pathway_genes <- NULL
  for (line in gmt) {
    split_line <- strsplit(line, "\t")[[1]]
    if (split_line[1] == pathway_name) {
      pathway_genes <- split_line[-c(1,2)]  # Assuming the first two columns are pathway name and description
      break
    }
  }
  
  # Check if pathway was found
  if (is.null(pathway_genes)) {
    stop("Pathway not found in GMT file.")
  }
  
  # Read survival data and sort patients by survival
  survival_df <- read.csv(survival_data, sep = "\t")
  sorted_patients <- survival_df[order(-survival_df$OS.from.enrollment..months.), 'Subject.ID']
  
  # Create a grid of feature plots
  timepoints <- c("N", "Y", "1", "2", "3", "4", "5", "R", "R2")
  # Initialize columns for mean expression at each timepoint in survival_df
  for(tp in timepoints) {
    survival_df[[paste0("Mean_Expr_", tp)]] <- NA_real_
  }
  
  for (patient in sorted_patients) {
    # print(patient)
    patient_specific_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$donor == patient,])
    if (length(patient_specific_cells) == 0) {
      next
    }
    patient_seurat_obj <- subset(seurat_obj, cells = patient_specific_cells)
    # Get the genes present in the Seurat object
    seurat_genes <- rownames(patient_seurat_obj@assays$RNA@counts)
    
    # Find intersection of pathway genes and Seurat object genes
    common_genes <- intersect(pathway_genes, seurat_genes)
    
    # Check if there are any common genes
    if (length(common_genes) == 0) {
      stop("None of the pathway genes are found in the Seurat object.")
    }
    
    # Calculate average expression of pathway genes per cell
    # Note: This assumes that the data is already normalized
    pathway_avg_expression <- apply(GetAssayData(patient_seurat_obj, assay = "RNA", slot = "data")[common_genes, ], 2, mean, na.rm = TRUE)
    
    # Check for cells at timepoint N
    timepoint_N_barcodes <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == "N", ])

    if (length(timepoint_N_barcodes) > 0) {
      # Calculate mean and standard deviation for cells at timepoint N
      mean_N <- mean(pathway_avg_expression[timepoint_N_barcodes])
      sd_N <- sd(pathway_avg_expression[timepoint_N_barcodes])

      # Scale based on timepoint N
      scaled_pathway_avg_expression <- (pathway_avg_expression - mean_N) / sd_N
    } else {
      # Use scale function if no cells at timepoint N
      scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    }

    pathway_avg_expression <- scaled_pathway_avg_expression
    names(pathway_avg_expression) <- colnames(patient_seurat_obj)
    
    # Add this as a metadata column
    patient_seurat_obj[["pathway_avg_expression"]] <- pathway_avg_expression
    
    for (tp in timepoints) {
      patient_tp_cells <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == tp, ])
      if (length(patient_tp_cells) > 0) {
        # Calculate mean expression for this timepoint
        mean_expression <- mean(pathway_avg_expression[patient_tp_cells], na.rm = TRUE)
      } else {
        mean_expression <- NA
      }
      # Update survival_df with mean expression for this patient and timepoint
      survival_df[survival_df$Subject.ID == patient, paste0("Mean_Expr_", tp)] <- mean_expression
    }
  }
  # Return the updated survival dataframe
  return(survival_df)
}

get_cluster_proportion <- function(seurat_metadata, patient_col, timepoint_col, cluster_col, clusters_of_interest, timepoint_of_interest) {
  # Filter metadata for the cluster of interest and the specified timepoint
  cluster_timepoint_data <- seurat_metadata[seurat_metadata[[cluster_col]] %in% clusters_of_interest & seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  timepoint_data <- seurat_metadata[seurat_metadata[[timepoint_col]] == timepoint_of_interest, ]
  # Group the filtered data by patient
  cluster_timpoint_grouped_data <- split(cluster_timepoint_data, cluster_timepoint_data[[patient_col]])
  timpoint_grouped_data <- split(timepoint_data, timepoint_data[[patient_col]])
  
  # Initialize lists to store Subject.ID and Cluster_Proportion_X
  subject_ids <- c()
  cluster_proportions <- c()
  
  # Loop through each patient
  for (key in names(cluster_timpoint_grouped_data)) {
    patient_cluster_timepoint_data <- cluster_timpoint_grouped_data[[key]]
    patient_timepoint_data <- timpoint_grouped_data[[key]]
    
    # Calculate the total number of cells for the patient
    patient_cluster_timepoint_cells <- nrow(patient_cluster_timepoint_data)
    patient_timepoint_cells <- nrow(patient_timepoint_data)
    
    # Add the total number of cells to the results if there are any cells
    if (patient_cluster_timepoint_cells > 0) {
      # Calculate the cluster proportion
      cluster_proportion <- patient_cluster_timepoint_cells / patient_timepoint_cells
      
      # Append the Subject.ID and Cluster_Proportion_X to the lists
      subject_ids <- c(subject_ids, key)
      cluster_proportions <- c(cluster_proportions, cluster_proportion)
    }
  }
  
  custom_column_name <- paste("Cluster_Proportion", timepoint_of_interest, sep = "_")
  # Create a data frame with Subject.ID and custom Cluster_Proportion_X column
  result_df <- data.frame(Subject.ID = subject_ids)
  result_df[[custom_column_name]] <- cluster_proportions
  
  return(result_df)
}

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/2TT_seurat_object_ITT.rds")
seurat_metadata <- seurat_object_all_cells@meta.data
# T1IFN resection type comparison
cluster_name <- "DC"
cluster_list <- c(29, 20)
# cluster_name <- "All_Cells"
# cluster_list <- unique(seurat_object_all_cells@meta.data$seurat_clusters)
seurat_object_subset <- subset(seurat_object_all_cells, subset = seurat_clusters %in% cluster_list)
seurat_object_subset <- NormalizeData(seurat_object_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"
pathway <- "GO_RESPONSE_TO_TYPE_I_INTERFERON"
# pathway <- "GO_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS"
t1ifn_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_subset, survival_data)

cox_input_df <- t1ifn_survival_df

timepoint_a <- "N"
timepoint_b <- "Y"


timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])


timepoint_b <- "1"

timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])


timepoint_b <- "2"

timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])


cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"  
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value
values_to_replace <- c("Biopsy")
new_value <- "Biopsy Only"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value
names(cox_input_df)[names(cox_input_df) == "Extension.of.resection..Biopsy.only..Partial..GTR."] <- "resection_type"

# Convert the DataFrame to long format
df_long <- pivot_longer(cox_input_df, cols = starts_with("Blood"), names_to = "count_type", values_to = "count_value")

df_long$count_type <- factor(df_long$count_type, levels = c("Blood_N_vs_Blood_Y_T1IFN_expression", "Blood_N_vs_Blood_1_T1IFN_expression", "Blood_N_vs_Blood_2_T1IFN_expression"))
# df_long$count_type <- factor(df_long$count_type, levels = c("Blood_N_vs_Blood_Y_IFL_expression", "Blood_N_vs_Blood_1_IFL_expression", "Blood_N_vs_Blood_2_IFL_expression"))
df_long$resection_type <- factor(df_long$resection_type, levels = c("Maximal Resection", "Biopsy Only"))

# Define custom colors
custom_colors <- c("Biopsy Only" = "red", "Maximal Resection" = "blue")

# Calculate p-values
results <- df_long %>%
  group_by(count_type) %>%
  summarise(
    p_value = {
      # Filtering the current data for necessary comparison
      current_data <- cur_data()
      if (all(c("Maximal Resection", "Biopsy Only") %in% unique(current_data$resection_type)) &&
          sum(current_data$resection_type == "Maximal Resection") >= 2 &&
          sum(current_data$resection_type == "Biopsy Only") >= 2) {
        # Perform Wilcoxon test
        test_result <- wilcox.test(count_value ~ resection_type, data = current_data)
        test_result$p.value
      } else {
        NA_real_  # Use NA_real_ for numerical context missing values
      }
    },
    .groups = 'drop'
  )

# Merge results back to the original data frame for plotting
df_long <- df_long %>%
  left_join(results, by = "count_type")

# Compute a position for the annotation just above the highest count_value per count_type
df_long <- df_long %>%
  group_by(count_type) %>%
  mutate(max_count_value = max(count_value, na.rm = TRUE)) %>%
  ungroup()

p <- ggplot(df_long, aes(x = count_type, y = count_value, fill = resection_type)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
  geom_jitter(position = position_dodge(width = 1.2)) +  # Correct position dodging for jitter
  facet_wrap(~count_type, scales = "free") +  # Create a separate plot for each count type
  labs(title = "IFL Signal Change by Resection Type", x = "Timepoint Comparison", y = "Signal Change") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        axis.title.x = element_text(size = 20),  # Adjust x axis label size
        axis.title.y = element_text(size = 20),  # Adjust y axis label size
        plot.title = element_text(size = 9),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),  # Add axis lines
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks.y = element_line()) +  # Add Y axis ticks
  scale_fill_manual(values = custom_colors) +
  geom_text(aes(y = Inf, label = ifelse(!is.na(p_value), sprintf("p = %.3f", p_value), "")), vjust = 1, color = "black", size = 7)

# # Print the plot
print(p)

plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/"
file_name <- paste0(cluster_name, "_T1IFN_with_cluster_proportion_Resection_Type_Comparison.pdf")
# file_name <- paste0(cluster_name, "_IFL_with_cluster_proportion_Resection_Type_Comparison.pdf")
ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 14, height = 4, dpi = 300, limitsize = FALSE)