library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)

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

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/2TT_T_Cells_seurat_object_ITT.rds")

include_cluster_proportion <- FALSE
patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathways <- c("GO_T_CELL_ACTIVATION", "GO_ADAPTIVE_IMMUNE_RESPONSE")
t_cell_subtypes <- c("all", "CD4", "CD8")

for (pathway in pathways) {
  for (t_cell_subtype in t_cell_subtypes) {
    if (t_cell_subtype == "all") {
      seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
      #######################################################################
      tcell_barcodes <- rownames(seurat_object_t_cell_subset@meta.data)
      temp <- seurat_object_all_cells@meta.data[tcell_barcodes,]
      tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
      #######################################################################
      print(tcell_cluster_list)
    } else if (t_cell_subtype == "CD4") {
      cd4_clusters <- c(1, 3, 6, 7, 9)
      seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cd4_clusters)
      seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
      tcell_cluster_list <- c(1,14)
    } else if (t_cell_subtype == "CD8") {
      cd8_clusters <- c(0, 5, 8, 10, 11, 14)
      seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cd8_clusters)
      seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
      tcell_cluster_list <- c(0,22)
    }
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
    cox_input_df <- updated_survival_df
    
    timepoint_a <- "N"
    timepoint_b <- "1"
    
    if (include_cluster_proportion) {
      timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
      timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
      cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
      cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
      cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
    } else {
      cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])
    }
    
    timepoint_b <- "2"
    
    if (include_cluster_proportion) {
      timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
      cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
      cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
    } else {
      cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] - cox_input_df[, paste0("Mean_Expr_", timepoint_a)])
    }
    
    names(cox_input_df)[names(cox_input_df) == "Extension.of.resection..Biopsy.only..Partial..GTR."] <- "resection_type"
    values_to_replace <- c("GTR", "Partial")
    new_value <- "Maximal Resection"  
    cox_input_df[, "resection_type"][cox_input_df[, "resection_type"] %in% values_to_replace] <- new_value
    values_to_replace <- c("Biopsy")
    new_value <- "Biopsy Only"
    cox_input_df[, "resection_type"][cox_input_df[, "resection_type"] %in% values_to_replace] <- new_value
    
    # Convert the DataFrame to long format
    df_long <- pivot_longer(cox_input_df, cols = starts_with("Blood"), names_to = "count_type", values_to = "count_value")
    
    df_long$count_type <- factor(df_long$count_type, levels = c("Blood_N_vs_Blood_1", "Blood_N_vs_Blood_2"))
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
      labs(title = paste0(pathway, " Signal Change by Resection Type in ", t_cell_subtype, " T Cells"), x = "Timepoint Comparison", y = "Signal Change") +
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
    if (include_cluster_proportion) {
      file_name <- paste0(pathway, "_", t_cell_subtype, "_T_Cells_with_cluster_proportion_Resection_Type_Comparison.pdf")
    } else {
      file_name <- paste0(pathway, "_", t_cell_subtype, "_T_Cells_without_cluster_proportion_Resection_Type_Comparison.pdf")
    }
    ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 14, height = 4, dpi = 300, limitsize = FALSE)
  }
}
######################################################################################################################################################
# instead of plotting change in expression, plot expression
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/2TT_seurat_object_ITT.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/2TT_T_Cells_seurat_object_ITT.rds")

include_cluster_proportion <- FALSE
patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathways <- c("GO_T_CELL_ACTIVATION", "GO_ADAPTIVE_IMMUNE_RESPONSE")
t_cell_subtypes <- c("all", "CD4", "CD8")

create_combined_plot <- function(df) {
  print(paste0("number of patients: ", nrow(df)))
  # Step 1: Filter columns related to Mean_Expr_
  expr_columns <- grep("Mean_Expr_", names(df), value = TRUE)
  
  # Step 2: Reshape the data to long format
  long_df <- df %>%
    pivot_longer(cols = expr_columns, names_to = "timepoint", values_to = "mean_expr") %>%
    mutate(timepoint = sub("Mean_Expr_", "", timepoint)) %>%
    filter(!timepoint %in% c("N"))  # Exclude timepoint "N"
  # filter(!timepoint %in% c("N", "3", "4", "5", "R","R2"))  # Exclude timepoint "N"
  
  # Calculate mean expression for each group at each timepoint
  summarized_df <- long_df %>%
    group_by(resection_type, timepoint) %>%
    summarize(mean_expr = median(mean_expr, na.rm = TRUE), .groups = 'drop') %>%
    filter(!is.nan(mean_expr))  # Ensure mean expressions are not NaN
  
  # Filtering to ensure we only include timepoints where both groups have data
  valid_timepoints <- summarized_df %>%
    group_by(timepoint) %>%
    filter(n() == 2) %>%
    pull(timepoint) %>%
    unique()
  
  summarized_df <- summarized_df %>%
    filter(timepoint %in% valid_timepoints)
  
  # Prepare data for plotting and statistical testing
  plot_df <- summarized_df %>%
    pivot_wider(names_from = resection_type, values_from = mean_expr)
  
  # Check for complete cases for paired t-test
  plot_df <- plot_df %>%
    filter(complete.cases(.))
  
  # Define custom colors
  custom_colors <- c("Biopsy Only" = "red", "Maximal Resection" = "blue")
  
  data = plot_df %>% gather(key = "resection_type", value = "mean_expr", -timepoint)
  
  data$resection_type <- factor(data$resection_type, levels = c("Maximal Resection", "Biopsy Only"))
  data$timepoint <- factor(data$timepoint, levels = c("Y", "1", "2", "3", "4", "5", "R", "R2"))
  
  # Statistical comparison using paired t-test
  if (nrow(plot_df) > 0) {
    t_test_result <- t.test(plot_df$`Biopsy Only`, plot_df$`Maximal Resection`, paired = TRUE)
    p_value <- t_test_result$p.value
    cat("P-value from paired t-test:", p_value, "\n")
  } else {
    p_value <- NA
    cat("Not enough data to perform a paired t-test.\n")
  }
  
  # Plotting
  p <- ggplot(data, aes(x = resection_type, y = mean_expr)) +
    geom_boxplot(aes(fill = resection_type), width = 0.2, position = position_dodge(width = 0.4), outlier.shape = NA) +  # Box plots without outliers
    # geom_line(data = data, aes(color = timepoint, group = timepoint), position = position_dodge(width = 0.75), size = 1) + # Adding lines between dots
    geom_line(data = data, aes(color = timepoint, group = timepoint), size = 1) + # Adding lines between dots
    # geom_point(data = data, aes(color = timepoint), position = position_dodge(width = 0.75), size = 5, alpha = 0.9) +  # Adding colored points
    geom_point(data = data, aes(color = timepoint), size = 3, alpha = 0.9) +  # Adding colored points
    scale_fill_manual(values = custom_colors) +
    scale_color_brewer(palette = "Set1", name = "Timepoint") +  # Colors for the dots
    labs(title = paste0("Comparison of Median Expression Between Groups", "\n", sprintf("p = %.9f", p_value)),
         x = "",
         y = "Median Expression per Timepoint",
         fill = "Resection Type") +
    theme_minimal() +
    theme(axis.line = element_line(color = "black"),
          axis.title.x = element_text(size = 10),  # Adjust x axis label size
          axis.title.y = element_text(size = 10),  # Adjust y axis label size
          plot.title = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),  # Add axis lines
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.ticks.y = element_line()) # Add Y axis ticks
  
  
  # Print the plot
  # print(p)
  return(p)
}


for (pathway in pathways) {
  for (t_cell_subtype in t_cell_subtypes) {
    if (t_cell_subtype == "all") {
      seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
      #######################################################################
      tcell_barcodes <- rownames(seurat_object_t_cell_subset@meta.data)
      temp <- seurat_object_all_cells@meta.data[tcell_barcodes,]
      tcell_cluster_list <- na.omit(as.numeric(as.character(unique(temp$seurat_clusters))))
      #######################################################################
      print(tcell_cluster_list)
    } else if (t_cell_subtype == "CD4") {
      cd4_clusters <- c(1, 3, 6, 7, 9)
      seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cd4_clusters)
      seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
      tcell_cluster_list <- c(1,14)
    } else if (t_cell_subtype == "CD8") {
      cd8_clusters <- c(0, 5, 8, 10, 11, 14)
      seurat_object_t_cell_subset <- subset(seurat_object_t_cells, subset = seurat_clusters %in% cd8_clusters)
      seurat_object_t_cell_subset <- NormalizeData(seurat_object_t_cell_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")
      tcell_cluster_list <- c(0,22)
    }
    
    updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cell_subset, survival_data)
    
    if (include_cluster_proportion) {
      # Create a grid of feature plots
      timepoints <- c("N", "Y", "1", "2", "3", "4", "5", "R", "R2")
      for (timepoint in timepoints) {
        timepoint_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint)
        updated_survival_df <- merge(updated_survival_df, timepoint_cluster_proportion_df, by = "Subject.ID", all.x = TRUE)
      }
      
      # List Mean_Exp and Cluster_Prop columns
      mean_exp_cols <- names(updated_survival_df)[grepl("Mean_Expr", names(updated_survival_df))]
      cluster_prop_cols <- names(updated_survival_df)[grepl("Cluster_Proportion", names(updated_survival_df))]
      
      # Ensure the number of Mean_Exp and Cluster_Prop columns match
      if (length(mean_exp_cols) == length(cluster_prop_cols)) {
        # Perform the multiplication
        for (i in seq_along(mean_exp_cols)) {
          updated_survival_df[[mean_exp_cols[i]]] <- updated_survival_df[[mean_exp_cols[i]]] * updated_survival_df[[cluster_prop_cols[i]]]
        }
      } else {
        stop("The number of Mean_Exp and Cluster_Prop columns does not match.")
      }
    }
    
    names(updated_survival_df)[names(updated_survival_df) == "Extension.of.resection..Biopsy.only..Partial..GTR."] <- "resection_type"
    values_to_replace <- c("GTR", "Partial")
    new_value <- "Maximal Resection"  
    updated_survival_df[, "resection_type"][updated_survival_df[, "resection_type"] %in% values_to_replace] <- new_value
    values_to_replace <- c("Biopsy")
    new_value <- "Biopsy Only"
    updated_survival_df[, "resection_type"][updated_survival_df[, "resection_type"] %in% values_to_replace] <- new_value
    
    # Convert the DataFrame to long format
    df_long <- pivot_longer(updated_survival_df, cols = starts_with("Mean"), names_to = "count_type", values_to = "count_value")
    
    # df_long$count_type <- factor(df_long$count_type, levels = c("Blood_N_vs_Blood_1", "Blood_N_vs_Blood_2"))
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
          unique_types <- unique(current_data$resection_type)
          
          # Check if exactly two unique resection types are present
          if (length(unique_types) == 2 && all(c("Maximal Resection", "Biopsy Only") %in% unique_types)) {
            # Ensure each group has at least two non-NA observations
            max_res_samples <- current_data$count_value[current_data$resection_type == "Maximal Resection"]
            biopsy_only_samples <- current_data$count_value[current_data$resection_type == "Biopsy Only"]
            
            if (sum(!is.na(max_res_samples)) >= 2 && sum(!is.na(biopsy_only_samples)) >= 2) {
              # Perform Wilcoxon test
              test_result <- wilcox.test(count_value ~ resection_type, data = current_data)
              if (!is.na(test_result$p.value)) {
                test_result$p.value
              } else {
                NA_real_
              }
            } else {
              NA_real_
            }
          } else {
            NA_real_
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
    
    p_each_timepoint <- ggplot(df_long, aes(x = count_type, y = count_value, fill = resection_type)) +
      geom_boxplot(width = 0.5, position = position_dodge(width = 1.2)) +  # Adjust dodge width for spacing
      geom_jitter(position = position_dodge(width = 1.2)) +  # Correct position dodging for jitter
      facet_wrap(~count_type, scales = "free") +  # Create a separate plot for each count type
      labs(title = paste0(pathway, " Signal Change by Resection Type in ", t_cell_subtype, " T Cells"), x = "Timepoint Comparison", y = "Signal Change") +
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
      scale_fill_manual(values = custom_colors)
    # geom_text(aes(y = Inf, label = ifelse(!is.na(p_value), sprintf("p = %.3f", p_value), "")), vjust = 1, color = "black", size = 7)
    
    # Print the plot
    # print(p_each_timepoint)
    
    p_combined <- create_combined_plot(updated_survival_df)
    
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/"
    if (include_cluster_proportion) {
      file_name_1 <- paste0(pathway, "_", t_cell_subtype, "_T_Cells_with_cluster_proportion_Resection_Type_Comparison_per_timepoint.pdf")
      file_name_2 <- paste0(pathway, "_", t_cell_subtype, "_T_Cells_with_cluster_proportion_Resection_Type_Comparison_combined_timepoints.pdf")
    } else {
      file_name_1 <- paste0(pathway, "_", t_cell_subtype, "_T_Cells_without_cluster_proportion_Resection_Type_Comparison_per_timepoint.pdf")
      file_name_2 <- paste0(pathway, "_", t_cell_subtype, "_T_Cells_without_cluster_proportion_Resection_Type_Comparison_combined_timepoints.pdf")
    }
    ggsave(paste0(plot_dir, file_name_1), p_each_timepoint, device = "pdf", width = 14, height = 10, dpi = 300, limitsize = FALSE)
    ggsave(paste0(plot_dir, file_name_2), p_combined, device = "pdf", width = 6, height = 6, dpi = 300, limitsize = FALSE)
  }
}