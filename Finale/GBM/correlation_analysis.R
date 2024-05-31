library(ggplot2)
library(survival)
library(Seurat)
library(ggfortify)
library("Rcpp")

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


# correlation with clonal replacement ratio
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data
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

t1ifn_timepoints <- c("N", "Y", "1", "2")
clonal_rep_timepoints <- c("Pre-TTF", "Post-TTF", "1", "2")
t1ifn_combinations <- apply(combn(t1ifn_timepoints, 2), 2, function(x) paste(x, collapse = ","))
clonal_rep_combinations <- apply(combn(clonal_rep_timepoints, 2), 2, function(x) paste(x, collapse = ","))

t1ifn_pathway <- "GO_RESPONSE_TO_TYPE_I_INTERFERON"
t1ifn_survival_df <- create_survival_data(gmt_file, t1ifn_pathway, seurat_object_subset, survival_data)
t1ifn_survival_df <- t1ifn_survival_df[t1ifn_survival_df$IDH.1.mutation == "negative",]

for (t1ifn_combination in t1ifn_combinations){
  for (clonal_rep_combination in clonal_rep_combinations){
    t1ifn_timepoint_a <- strsplit(t1ifn_combination, split = ",")[[1]][1]
    t1ifn_timepoint_b <- strsplit(t1ifn_combination, split = ",")[[1]][2]
    
    clonal_rep_timepoint_a <- strsplit(clonal_rep_combination, split = ",")[[1]][1]
    clonal_rep_timepoint_b <- strsplit(clonal_rep_combination, split = ",")[[1]][2]
    
    cox_input_df <- t1ifn_survival_df
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, t1ifn_timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, t1ifn_timepoint_b)
    cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_b)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_b)]), ]
    # cox_input_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_a)])
    cox_input_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_a)])
    
    # read cloanl replacement ratio file
    clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
    clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", clonal_rep_timepoint_a, "_sc_vs_Blood_", clonal_rep_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
    
    merged_df <- merge(cox_input_df[, c("Subject.ID", "OS.from.enrollment..months.", paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression"))], clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    # correlation_pearson <- cor(merged_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")], merged_df$top_10, use = "complete.obs")
    correlation_pearson <-cor.test(merged_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")], merged_df$top_10, method = "pearson")
    print(correlation_pearson)
    
    x = paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")
    y = "top_10"
    # Create a scatter plot with ggplot2
    p <- ggplot(merged_df, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(colour = "black") +  # Adding black outline to points
      geom_smooth(method = "lm", formula = y ~ x, col = "red") +  # Linear regression line
      # ylim(-0.005, 0.037) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  # Black outline
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank()
      ) +
      labs(
        title = sprintf("Scatter Plot with Pearson r = %.2f, p-value = %.6f",
                        correlation_pearson$estimate, correlation_pearson$p.value),
        x = paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression"),  # Removing x-axis label
        y = paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Replacement")   # Removing y-axis label
      )
    
    # Print the plot
    print(p)
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/T1IFN_Expression_vs_Clonal_Replacement/"
    x_axis_lab <- paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")
    y_axis_lab <- paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Replacement")
    file_name <- paste0(cluster_name, "_", x_axis_lab, "_vs_", y_axis_lab, "_without_cluster_proportion.png")
    ggsave(paste0(plot_dir, file_name), p, width = 6, height = 4, dpi = 300)
  }
}



###########################################################################################################################################################################################
# correlation between T CELL ACTIVATION expression and T CELL Replacement

# correlation with clonal replacement ratio
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_T_Cells_seurat_object_GBM.rds")
seurat_object_t_cells <- NormalizeData(seurat_object_t_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

tcell_timepoints <- c("N", "Y", "1", "2")
clonal_rep_timepoints <- c("Pre-TTF", "Post-TTF", "1", "2")
tcell_combinations <- apply(combn(tcell_timepoints, 2), 2, function(x) paste(x, collapse = ","))
clonal_rep_combinations <- apply(combn(clonal_rep_timepoints, 2), 2, function(x) paste(x, collapse = ","))

pathway <- "GO_T_CELL_ACTIVATION"
tcell_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cells, survival_data)
tcell_survival_df <- tcell_survival_df[tcell_survival_df$IDH.1.mutation == "negative",]

#######################################################################
cluster_name <- "T_Cells"
tcell_barcodes <- rownames(seurat_object_t_cells@meta.data)
temp <- seurat_object_all_cells@meta.data[tcell_barcodes,]
tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
#######################################################################

for (tcell_combination in tcell_combinations){
  for (clonal_rep_combination in clonal_rep_combinations){
    tcell_timepoint_a <- strsplit(tcell_combination, split = ",")[[1]][1]
    tcell_timepoint_b <- strsplit(tcell_combination, split = ",")[[1]][2]
    
    clonal_rep_timepoint_a <- strsplit(clonal_rep_combination, split = ",")[[1]][1]
    clonal_rep_timepoint_b <- strsplit(clonal_rep_combination, split = ",")[[1]][2]
    
    cox_input_df <- tcell_survival_df
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, tcell_timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, tcell_timepoint_b)
    cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_b)]), ]
    cox_input_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")] <- (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_a)])
    # cox_input_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")] <- (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)])
    
    # read cloanl replacement ratio file
    clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
    clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", clonal_rep_timepoint_a, "_sc_vs_Blood_", clonal_rep_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
    
    merged_df <- merge(cox_input_df[, c("Subject.ID", "OS.from.enrollment..months.", paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression"))], clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    correlation_pearson <-cor.test(merged_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")], merged_df$top_10, method = "pearson")
    print(correlation_pearson)
    
    x = paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")
    y = "top_10"
    # Create a scatter plot with ggplot2
    p <- ggplot(merged_df, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(colour = "black") +  # Adding black outline to points
      geom_smooth(method = "lm", formula = y ~ x, col = "red") +  # Linear regression line
      # ylim(-0.005, 0.037) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  # Black outline
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank()
      ) +
      labs(
        title = sprintf("Scatter Plot with Pearson r = %.2f, p-value = %.6f",
                        correlation_pearson$estimate, correlation_pearson$p.value),
        x = paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression"),  # Removing x-axis label
        y = paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Replacement")   # Removing y-axis label
      )
    
    # Print the plot
    print(p)
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/TCA_Expression_vs_Clonal_Replacement/"
    x_axis_lab <- paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")
    y_axis_lab <- paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Replacement")
    file_name <- paste0(cluster_name, "_", x_axis_lab, "_vs_", y_axis_lab, "_with_cluster_proportion.png")
    ggsave(paste0(plot_dir, file_name), p, width = 6, height = 4, dpi = 300)
  }
}


############################################################################################################
# correlation between DC-T1IFN to the N vs Y clonal expansion


seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data
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

t1ifn_timepoints <- c("N", "Y", "1", "2")
clonal_rep_timepoints <- c("Pre-TTF", "Post-TTF")
t1ifn_combinations <- apply(combn(t1ifn_timepoints, 2), 2, function(x) paste(x, collapse = ","))
clonal_rep_combinations <- apply(combn(clonal_rep_timepoints, 2), 2, function(x) paste(x, collapse = ","))

t1ifn_pathway <- "GO_RESPONSE_TO_TYPE_I_INTERFERON"
t1ifn_survival_df <- create_survival_data(gmt_file, t1ifn_pathway, seurat_object_subset, survival_data)
t1ifn_survival_df <- t1ifn_survival_df[t1ifn_survival_df$IDH.1.mutation == "negative",]

for (t1ifn_combination in t1ifn_combinations){
  for (clonal_rep_combination in clonal_rep_combinations){
    t1ifn_timepoint_a <- strsplit(t1ifn_combination, split = ",")[[1]][1]
    t1ifn_timepoint_b <- strsplit(t1ifn_combination, split = ",")[[1]][2]
    
    clonal_rep_timepoint_a <- strsplit(clonal_rep_combination, split = ",")[[1]][1]
    clonal_rep_timepoint_b <- strsplit(clonal_rep_combination, split = ",")[[1]][2]
    
    cox_input_df <- t1ifn_survival_df
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, t1ifn_timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, t1ifn_timepoint_b)
    cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_b)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_b)]), ]
    # cox_input_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", t1ifn_timepoint_a)])
    cox_input_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", t1ifn_timepoint_a)])
    
    # read cloanl replacement ratio file
    clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
    clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", clonal_rep_timepoint_a, "_sc_vs_Blood_", clonal_rep_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
    
    merged_df <- merge(cox_input_df[, c("Subject.ID", "OS.from.enrollment..months.", paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression"))], clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    # correlation_pearson <- cor(merged_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")], merged_df$top_10, use = "complete.obs")
    correlation_pearson <-cor.test(merged_df[, paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")], merged_df$clonal_expansion, method = "pearson")
    print(correlation_pearson)
    
    x = paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")
    y = "clonal_expansion"
    # Create a scatter plot with ggplot2
    p <- ggplot(merged_df, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(colour = "black") +  # Adding black outline to points
      geom_smooth(method = "lm", formula = y ~ x, col = "red") +  # Linear regression line
      # ylim(-0.005, 0.037) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  # Black outline
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank()
      ) +
      labs(
        title = sprintf("Scatter Plot with Pearson r = %.2f, p-value = %.6f",
                        correlation_pearson$estimate, correlation_pearson$p.value),
        x = paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression"),  # Removing x-axis label
        y = paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Diversity")   # Removing y-axis label
      )
    
    # Print the plot
    print(p)
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/T1IFN_Expression_vs_Clonal_Diversity/"
    x_axis_lab <- paste0("Blood_", t1ifn_timepoint_a, "_vs_Blood_", t1ifn_timepoint_b, "_T1IFN_expression")
    y_axis_lab <- paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Diversity")
    file_name <- paste0(cluster_name, "_", x_axis_lab, "_vs_", y_axis_lab, "_without_cluster_proportion.png")
    ggsave(paste0(plot_dir, file_name), p, width = 6, height = 4, dpi = 300)
  }
}


###########################################################################################################################################################################################
# correlation between T CELL ACTIVATION expression and T CELL Diversity
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_T_Cells_seurat_object_GBM.rds")
seurat_object_t_cells <- NormalizeData(seurat_object_t_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

tcell_timepoints <- c("N", "Y", "1", "2")
clonal_rep_timepoints <- c("Pre-TTF", "Post-TTF")
tcell_combinations <- apply(combn(tcell_timepoints, 2), 2, function(x) paste(x, collapse = ","))
clonal_rep_combinations <- apply(combn(clonal_rep_timepoints, 2), 2, function(x) paste(x, collapse = ","))

pathway <- "GO_T_CELL_ACTIVATION"
tcell_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cells, survival_data)
tcell_survival_df <- tcell_survival_df[tcell_survival_df$IDH.1.mutation == "negative",]

#######################################################################
cluster_name <- "T_Cells"
tcell_barcodes <- rownames(seurat_object_t_cells@meta.data)
temp <- seurat_object_all_cells@meta.data[tcell_barcodes,]
tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
#######################################################################

for (tcell_combination in tcell_combinations){
  for (clonal_rep_combination in clonal_rep_combinations){
    tcell_timepoint_a <- strsplit(tcell_combination, split = ",")[[1]][1]
    tcell_timepoint_b <- strsplit(tcell_combination, split = ",")[[1]][2]
    
    clonal_rep_timepoint_a <- strsplit(clonal_rep_combination, split = ",")[[1]][1]
    clonal_rep_timepoint_b <- strsplit(clonal_rep_combination, split = ",")[[1]][2]
    
    cox_input_df <- tcell_survival_df
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, tcell_timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, tcell_timepoint_b)
    cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_b)]), ]
    cox_input_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")] <- (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_a)])
    # cox_input_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")] <- (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)])
    
    # read cloanl replacement ratio file
    clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
    clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", clonal_rep_timepoint_a, "_sc_vs_Blood_", clonal_rep_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
    
    merged_df <- merge(cox_input_df[, c("Subject.ID", "OS.from.enrollment..months.", paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression"))], clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    correlation_pearson <-cor.test(merged_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")], merged_df$clonal_expansion, method = "pearson")
    print(correlation_pearson)
    
    x = paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")
    y = "clonal_expansion"
    # Create a scatter plot with ggplot2
    p <- ggplot(merged_df, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(colour = "black") +  # Adding black outline to points
      geom_smooth(method = "lm", formula = y ~ x, col = "red") +  # Linear regression line
      # ylim(-0.005, 0.037) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  # Black outline
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank()
      ) +
      labs(
        title = sprintf("Scatter Plot with Pearson r = %.2f, p-value = %.6f",
                        correlation_pearson$estimate, correlation_pearson$p.value),
        x = paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression"),  # Removing x-axis label
        y = paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Diversity")   # Removing y-axis label
      )
    
    # Print the plot
    print(p)
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/TCA_Expression_vs_Clonal_Diversity/"
    x_axis_lab <- paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_TCA_expression")
    y_axis_lab <- paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Diversity")
    file_name <- paste0(cluster_name, "_", x_axis_lab, "_vs_", y_axis_lab, "_with_cluster_proportion.png")
    ggsave(paste0(plot_dir, file_name), p, width = 6, height = 4, dpi = 300)
  }
}


###########################################################################################################################################################################################
# correlation between ADAPTIVE IMMUNE RESPOSNSE expression and T CELL Diversity
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_metadata <- seurat_object_all_cells@meta.data

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_T_Cells_seurat_object_GBM.rds")
seurat_object_t_cells <- NormalizeData(seurat_object_t_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

tcell_timepoints <- c("N", "Y", "1", "2")
clonal_rep_timepoints <- c("Pre-TTF", "Post-TTF")
tcell_combinations <- apply(combn(tcell_timepoints, 2), 2, function(x) paste(x, collapse = ","))
clonal_rep_combinations <- apply(combn(clonal_rep_timepoints, 2), 2, function(x) paste(x, collapse = ","))

pathway <- "GO_ADAPTIVE_IMMUNE_RESPONSE"
tcell_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cells, survival_data)
tcell_survival_df <- tcell_survival_df[tcell_survival_df$IDH.1.mutation == "negative",]

#######################################################################
cluster_name <- "T_Cells"
tcell_barcodes <- rownames(seurat_object_t_cells@meta.data)
temp <- seurat_object_all_cells@meta.data[tcell_barcodes,]
tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
#######################################################################

for (tcell_combination in tcell_combinations){
  for (clonal_rep_combination in clonal_rep_combinations){
    tcell_timepoint_a <- strsplit(tcell_combination, split = ",")[[1]][1]
    tcell_timepoint_b <- strsplit(tcell_combination, split = ",")[[1]][2]
    
    clonal_rep_timepoint_a <- strsplit(clonal_rep_combination, split = ",")[[1]][1]
    clonal_rep_timepoint_b <- strsplit(clonal_rep_combination, split = ",")[[1]][2]
    
    cox_input_df <- tcell_survival_df
    timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, tcell_timepoint_a)
    timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, tcell_timepoint_b)
    cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_a)]), ]
    cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_b)]), ]
    cox_input_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression")] <- (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", tcell_timepoint_a)])
    # cox_input_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression")] <- (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", tcell_timepoint_a)])
    
    # read cloanl replacement ratio file
    clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
    clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", clonal_rep_timepoint_a, "_sc_vs_Blood_", clonal_rep_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
    
    merged_df <- merge(cox_input_df[, c("Subject.ID", "OS.from.enrollment..months.", paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression"))], clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    correlation_pearson <-cor.test(merged_df[, paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression")], merged_df$clonal_expansion, method = "pearson")
    print(correlation_pearson)
    
    x = paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression")
    y = "clonal_expansion"
    # Create a scatter plot with ggplot2
    p <- ggplot(merged_df, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(colour = "black") +  # Adding black outline to points
      geom_smooth(method = "lm", formula = y ~ x, col = "red") +  # Linear regression line
      # ylim(-0.005, 0.037) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  # Black outline
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank()
      ) +
      labs(
        title = sprintf("Scatter Plot with Pearson r = %.2f, p-value = %.6f",
                        correlation_pearson$estimate, correlation_pearson$p.value),
        x = paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression"),  # Removing x-axis label
        y = paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Diversity")   # Removing y-axis label
      )
    
    # Print the plot
    print(p)
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/AIR_Expression_vs_Clonal_Diversity/"
    x_axis_lab <- paste0("Blood_", tcell_timepoint_a, "_vs_Blood_", tcell_timepoint_b, "_AIR_expression")
    y_axis_lab <- paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Diversity")
    file_name <- paste0(cluster_name, "_", x_axis_lab, "_vs_", y_axis_lab, "_with_cluster_proportion.png")
    ggsave(paste0(plot_dir, file_name), p, width = 6, height = 4, dpi = 300)
  }
}


###########################################################################################################################################################################################
# correlation between T Cell Replacement expression and T Cell Diversity

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")
survival_df <- survival_df[survival_df$IDH.1.mutation == "negative",]

clonal_rep_timepoints <- c("Pre-TTF", "Post-TTF", "1", "2")
clonal_div_timepoints <- c("Pre-TTF", "Post-TTF")
clonal_rep_combinations <- apply(combn(clonal_rep_timepoints, 2), 2, function(x) paste(x, collapse = ","))
clonal_div_combinations <- apply(combn(clonal_div_timepoints, 2), 2, function(x) paste(x, collapse = ","))

for (clonal_rep_combination in clonal_rep_combinations){
  for (clonal_div_combination in clonal_div_combinations){
    clonal_rep_timepoint_a <- strsplit(clonal_rep_combination, split = ",")[[1]][1]
    clonal_rep_timepoint_b <- strsplit(clonal_rep_combination, split = ",")[[1]][2]
    
    clonal_div_timepoint_a <- strsplit(clonal_div_combination, split = ",")[[1]][1]
    clonal_div_timepoint_b <- strsplit(clonal_div_combination, split = ",")[[1]][2]
    
    cox_input_df <- survival_df
    
    # read cloanl replacement ratio file
    clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
    clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", clonal_rep_timepoint_a, "_sc_vs_Blood_", clonal_rep_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
    
    merged_df <- merge(cox_input_df[, c("Subject.ID", "OS.from.enrollment..months.")], clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    # read cloanl diversity ratio file
    clonal_diversity_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
    clonal_diversity_ratio <- read.table(paste0(clonal_diversity_dir, "Blood_", clonal_div_timepoint_a, "_sc_vs_Blood_", clonal_div_timepoint_b, "_sc.txt"), sep = "\t")
    rownames(clonal_diversity_ratio) <- as.numeric(sub("^p", "", rownames(clonal_diversity_ratio)))
    
    merged_df <- merge(merged_df, clonal_diversity_ratio, by.x = "Subject.ID", by.y = "row.names")
    
    correlation_pearson <-cor.test(merged_df$top_10, merged_df$clonal_expansion, method = "pearson")
    print(correlation_pearson)
    
    x = "top_10"
    y = "clonal_expansion"
    # Create a scatter plot with ggplot2
    p <- ggplot(merged_df, aes(x = !!sym(x), y = !!sym(y))) +
      geom_point(colour = "black") +  # Adding black outline to points
      geom_smooth(method = "lm", formula = y ~ x, col = "red") +  # Linear regression line
      # ylim(-0.005, 0.037) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  # Black outline
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank()
      ) +
      labs(
        title = sprintf("Scatter Plot with Pearson r = %.2f, p-value = %.6f",
                        correlation_pearson$estimate, correlation_pearson$p.value),
        x = paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Replacement"),  # Removing x-axis label
        y = paste0("Blood_", clonal_div_timepoint_a, "_vs_Blood_", clonal_div_timepoint_b, "_Clonal_Diversity")   # Removing y-axis label
      )
    
    # Print the plot
    print(p)
    plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/Clonal_Replacement_vs_Clonal_Diversity/"
    x_axis_lab <- paste0("Blood_", clonal_rep_timepoint_a, "_vs_Blood_", clonal_rep_timepoint_b, "_Clonal_Replacement")
    y_axis_lab <- paste0("Blood_", clonal_div_timepoint_a, "_vs_Blood_", clonal_div_timepoint_b, "_Clonal_Diversity")
    file_name <- paste0(x_axis_lab, "_vs_", y_axis_lab, ".png")
    ggsave(paste0(plot_dir, file_name), p, width = 6, height = 4, dpi = 300)
  }
}