library(ggplot2)
library(survival)
library(Seurat)
library(ggfortify)
library("Rcpp")

get_cox = function(cox_input_subset_df, time_column, event_column, covariates_list)
{
  df = cox_input_subset_df
  df =  df[, c(time_column, event_column, covariates_list)]
  # df[, new_covariate] <- scale(df[, new_covariate])
  
  res.cox <- NULL
  
  formula_string = paste0("Surv(", time_column, ", ", event_column, ") ~ ", paste(covariates_list, collapse = " + "))
  new_formula <- as.formula(formula_string)
  print(new_formula)
  tryCatch(
    {
      res.cox <- coxph(new_formula, data = df, control = coxph.control(iter.max=20))
    },
    error = function(e_outer) {
      cat("Outer Error: ", conditionMessage(e_outer), "\n")
    },
    warning = function(w_outer) {
      if (grepl("Ran out of iterations and did not converge", conditionMessage(w_outer))) {
        cat("Caught the specific warning: Ran out of iterations and did not converge\n")
        tryCatch(
          {
            res.cox <- coxph(new_formula, data = df, control = coxph.control(iter.max=1000))
          },
          error = function(e_inner) {
            cat("Inner Error: ", conditionMessage(e_inner), "\n")
          },
          warning = function(w_inner) {
            warning(w_inner)
          }
        )
      } else {
        # Handle other warnings if needed
        warning(w_outer)
      }
    }
  )
  
  return(list(cox = res.cox))
  
}

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")

timepoint_a <- "Post-TTF"
timepoint_b <- "2"
top_n <- "top_20"

# read cloanl replacement ratio file
in_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
clonal_replacement_ratio <- read.table(paste0(in_dir, paste0("Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt")), sep = "\t")
rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))

merged_df <- merge(survival_df, clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")

cox_input_df <- merged_df
cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value
cox_input_df <- cox_input_df[!is.na(cox_input_df[, top_n]), ]

# Define other covariates excluding 'top_10' which is always included
other_covariates <- c("Age", "Sex", "IDH.1.mutation", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")

concordance_list <- c()
hr_list <- c()
hr_upper_list <- c()
hr_lower_list <- c()
p_val_list <- c()
covariates_combination_list <- c()
expression_timepoint_list <- c()
cox_result_output_list <- list()

# Generate all combinations of the other covariates
for (i in 0:length(other_covariates)) {
  combinations <- combn(other_covariates, i, simplify = FALSE)
  for (combo in combinations) {
    # Always include 'top_10'
    covariates <- c(combo, top_n)
    # Run Cox analysis
    cox_ph_result <- get_cox(cox_input_df, "OS.from.enrollment..months.", "Dead", covariates)
    if (!is.null(cox_ph_result$cox)) {
      temp <- summary(cox_ph_result$cox)
      hr_list <- c(hr_list, temp$coefficients[top_n, "exp(coef)"])
      p_val_list <- c(p_val_list, temp$coefficients[top_n, "Pr(>|z|)"])
      hr_upper_list <- c(hr_upper_list, temp$conf.int[top_n, "upper .95"])
      hr_lower_list <- c(hr_lower_list, temp$conf.int[top_n, "lower .95"])
      concordance_list <- c(concordance_list, temp$concordance[["C"]])
      covariates_combination_list <- c(covariates_combination_list, paste(covariates, collapse = ", "))
    }
  }
}



df <- data.frame(covariate_combination = covariates_combination_list,
                 p_value = p_val_list,
                 hazard_ratio = hr_list,
                 hazard_ratio_upper_CI = hr_upper_list,
                 hazard_ratio_lower_CI = hr_lower_list,
                 concordance = concordance_list)


# only clinical covariates
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")

cox_input_df <- survival_df
cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value

covariates <- c("Age", "Sex", "IDH.1.mutation", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")
cox_ph_result <- get_cox(cox_input_df, "OS.from.enrollment..months.", "Dead", covariates)
summary(cox_ph_result$cox)


# calculating p-value for increase in corcordance

library(survival)
library(survminer)

# Function to calculate the p-value for survival analysis
calculate_p_value_survival <- function(time, status, covariates, new_covariate, n_permutations = 10000) {
  # Create the initial survival model with the new covariate
  original_data <- data.frame(time = time, status = status, covariates, new_covariate)
  formula <- as.formula(paste("Surv(time, status) ~", paste(names(covariates), collapse = "+"), "+ new_covariate"))
  original_model <- coxph(formula, data = original_data)
  orig_temp <- summary(original_model)
  original_concordance <- orig_temp$concordance[["C"]]
  
  # Array to store concordances from permutations
  permuted_concordances <- numeric(n_permutations)
  error_count <- 0  # Initialize error counter
  
  # Perform permutations
  set.seed(123)  # For reproducibility
  for (i in 1:n_permutations) {
    repeat {
      shuffled_covariate <- sample(original_data$new_covariate)
      if (!all(shuffled_covariate == original_data$new_covariate)) break
    }
    
    permuted_data <- original_data
    permuted_data$new_covariate <- shuffled_covariate
    
    # Attempt to fit model and calculate concordance
    permuted_result <- tryCatch({
      permuted_model <- coxph(formula, data = permuted_data)
      perm_temp <- summary(permuted_model)
      perm_temp$concordance[["C"]]
    }, error = function(e) {
      error_count <<- error_count + 1  # Increment error counter on failure
      NA  # Return NA on error
    })
    
    permuted_concordances[i] <- permuted_result
  }
  
  # Filter out NA values due to errors
  valid_concordances <- permuted_concordances[!is.na(permuted_concordances)]
  
  # Calculate p-value: proportion of valid permuted concordances >= actual concordance
  if (!is.na(original_concordance)) {
    p_value <- mean(valid_concordances >= original_concordance)
  } else {
    p_value <- NA  # If original concordance failed, p-value is undefined
  }
  
  # Return results including the error count
  list(actual_concordance = original_concordance, p_value = p_value, errors = error_count)
}

# Example usage:
# Assuming 'time', 'status', 'covariates', and 'new_covariate' are appropriately defined in your environment
# For instance:
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")

timepoint_a <- "Pre-TTF"
timepoint_b <- "2"
top_n <- "top_20"

# read cloanl replacement ratio file
in_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
clonal_replacement_ratio <- read.table(paste0(in_dir, paste0("Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt")), sep = "\t")
rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))

merged_df <- merge(survival_df, clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")

cox_input_df <- merged_df
cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value
cox_input_df <- cox_input_df[!is.na(cox_input_df[, top_n]), ]

time <- cox_input_df$OS.from.enrollment..months.
status <- cox_input_df$Dead
covariates <- cox_input_df[, c("Age", "Sex", "IDH.1.mutation", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")]
new_covariate <- cox_input_df[, c(top_n)]
result <- calculate_p_value_survival(time, status, covariates, new_covariate)

# Print result
print(result)
###############################################################################################################################################################################
# T1IFN Survival Analysis

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
    
    # # Check for cells at timepoint N
    # timepoint_N_barcodes <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == "N", ])
    # 
    # if (length(timepoint_N_barcodes) > 0) {
    #   # Calculate mean and standard deviation for cells at timepoint N
    #   mean_N <- mean(pathway_avg_expression[timepoint_N_barcodes])
    #   sd_N <- sd(pathway_avg_expression[timepoint_N_barcodes])
    #   
    #   # Scale based on timepoint N
    #   scaled_pathway_avg_expression <- (pathway_avg_expression - mean_N) / sd_N
    # } else {
    #   # Use scale function if no cells at timepoint N
    #   scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    # }
    # 
    # pathway_avg_expression <- scaled_pathway_avg_expression
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
cluster_name <- "DC"
cluster_list <- c(29, 20)
# cluster_name <- "Monocytes"
# cluster_list <- c(8, 11, 18, 2)
# cluster_name <- "NK"
# cluster_list <- c(3, 17)
seurat_object_subset <- subset(seurat_object_all_cells, subset = seurat_clusters %in% cluster_list)
seurat_object_subset <- NormalizeData(seurat_object_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"
pathway <- "GO_RESPONSE_TO_TYPE_I_INTERFERON"
t1ifn_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_subset, survival_data)

timepoint_a <- "N"
timepoint_b <- "1"

cox_input_df <- t1ifn_survival_df
timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", timepoint_a)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", timepoint_b)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]), ]
cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])

cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value

# Define other covariates excluding paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression") which is always included
other_covariates <- c("Age", "Sex", "IDH.1.mutation", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")

concordance_list <- c()
hr_list <- c()
hr_upper_list <- c()
hr_lower_list <- c()
p_val_list <- c()
covariates_combination_list <- c()
expression_timepoint_list <- c()
cox_result_output_list <- list()

constant_covariate <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")

# Generate all combinations of the other covariates
for (i in 0:length(other_covariates)) {
  combinations <- combn(other_covariates, i, simplify = FALSE)
  for (combo in combinations) {
    # Always include paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")
    covariates <- c(combo, constant_covariate)
    # Run Cox analysis
    cox_ph_result <- get_cox(cox_input_df, "OS.from.enrollment..months.", "Dead", covariates)
    if (!is.null(cox_ph_result$cox)) {
      temp <- summary(cox_ph_result$cox)
      hr_list <- c(hr_list, temp$coefficients[constant_covariate, "exp(coef)"])
      p_val_list <- c(p_val_list, temp$coefficients[constant_covariate, "Pr(>|z|)"])
      hr_upper_list <- c(hr_upper_list, temp$conf.int[constant_covariate, "upper .95"])
      hr_lower_list <- c(hr_lower_list, temp$conf.int[constant_covariate, "lower .95"])
      concordance_list <- c(concordance_list, temp$concordance[["C"]])
      covariates_combination_list <- c(covariates_combination_list, paste(covariates, collapse = ", "))
    }
  }
}



df <- data.frame(covariate_combination = covariates_combination_list,
                 p_value = p_val_list,
                 hazard_ratio = hr_list,
                 hazard_ratio_upper_CI = hr_upper_list,
                 hazard_ratio_lower_CI = hr_lower_list,
                 concordance = concordance_list)


##########################################################################################################################
cluster_name <- "DC"
cluster_list <- c(29, 20)
seurat_object_subset <- subset(seurat_object_all_cells, subset = seurat_clusters %in% cluster_list)
seurat_object_subset <- NormalizeData(seurat_object_subset, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"
pathway <- "GO_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS"
t1ifn_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_subset, survival_data)

timepoint_a <- "N"
timepoint_b <- "1"

cox_input_df <- t1ifn_survival_df
timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_a)
timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, cluster_list, timepoint_b)
cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", timepoint_a)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", timepoint_b)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]), ]
cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])

cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value

# Define other covariates excluding paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression") which is always included
other_covariates <- c("Age", "Sex", "IDH.1.mutation", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")

concordance_list <- c()
hr_list <- c()
hr_upper_list <- c()
hr_lower_list <- c()
p_val_list <- c()
covariates_combination_list <- c()
expression_timepoint_list <- c()
cox_result_output_list <- list()

constant_covariate <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_IFL_expression")

# Generate all combinations of the other covariates
for (i in 0:length(other_covariates)) {
  combinations <- combn(other_covariates, i, simplify = FALSE)
  for (combo in combinations) {
    # Always include paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")
    covariates <- c(combo, constant_covariate)
    # Run Cox analysis
    cox_ph_result <- get_cox(cox_input_df, "OS.from.enrollment..months.", "Dead", covariates)
    if (!is.null(cox_ph_result$cox)) {
      temp <- summary(cox_ph_result$cox)
      hr_list <- c(hr_list, temp$coefficients[constant_covariate, "exp(coef)"])
      p_val_list <- c(p_val_list, temp$coefficients[constant_covariate, "Pr(>|z|)"])
      hr_upper_list <- c(hr_upper_list, temp$conf.int[constant_covariate, "upper .95"])
      hr_lower_list <- c(hr_lower_list, temp$conf.int[constant_covariate, "lower .95"])
      concordance_list <- c(concordance_list, temp$concordance[["C"]])
      covariates_combination_list <- c(covariates_combination_list, paste(covariates, collapse = ", "))
    }
  }
}



df <- data.frame(covariate_combination = covariates_combination_list,
                 p_value = p_val_list,
                 hazard_ratio = hr_list,
                 hazard_ratio_upper_CI = hr_upper_list,
                 hazard_ratio_lower_CI = hr_lower_list,
                 concordance = concordance_list)


#######################################################################################################################################################################
# T cell Activation Survival Analysis

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
seurat_object_t_cells <- NormalizeData(seurat_object_t_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

patient_col <- "donor"
timepoint_col <- "TimePoint"
cluster_col <- "seurat_clusters"

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"

pathway <- "GO_T_CELL_ACTIVATION"

#######################################################################
cluster_name <- "T_Cells"
tcell_barcodes <- rownames(seurat_object_t_cells@meta.data)
temp <- seurat_object_all_cells@meta.data[tcell_barcodes,]
tcell_cluster_list <- as.numeric(as.character(unique(temp$seurat_clusters)))
#######################################################################

updated_survival_df <- create_survival_data(gmt_file, pathway, seurat_object_t_cells, survival_data) 

timepoint_a <- "N"
timepoint_b <- "2"

cox_input_df <- updated_survival_df
timepoint_a_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_a)
timepoint_b_cluster_proportion_df <- get_cluster_proportion(seurat_metadata, patient_col, timepoint_col, cluster_col, tcell_cluster_list, timepoint_b)
cox_input_df <- merge(cox_input_df, timepoint_a_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- merge(cox_input_df, timepoint_b_cluster_proportion_df, by = "Subject.ID")
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", timepoint_a)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Mean_Expr_", timepoint_b)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)]), ]
cox_input_df <- cox_input_df[!is.na(cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]), ]
cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_TCA_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)] * cox_input_df[, paste0("Cluster_Proportion_", timepoint_a)])
# cox_input_df[, paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_TCA_expression")] <- (cox_input_df[, paste0("Mean_Expr_", timepoint_b)]) - (cox_input_df[, paste0("Mean_Expr_", timepoint_a)])

cox_input_df$Dead <- ifelse(cox_input_df$Dead == "N", 0, ifelse(cox_input_df$Dead == "Y", 1, cox_input_df$Dead))
cox_input_df$Dead <- as.numeric(cox_input_df$Dead) 
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"
cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."][cox_input_df[, "Extension.of.resection..Biopsy.only..Partial..GTR."] %in% values_to_replace] <- new_value


# Define other covariates excluding 'clonal_expansion' which is always included
other_covariates <- c("Age", "Sex", "IDH.1.mutation", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")
constant_covariate <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_TCA_expression")

concordance_list <- c()
hr_list <- c()
hr_upper_list <- c()
hr_lower_list <- c()
p_val_list <- c()
covariates_combination_list <- c()
expression_timepoint_list <- c()
cox_result_output_list <- list()

# Generate all combinations of the other covariates
for (i in 0:length(other_covariates)) {
  combinations <- combn(other_covariates, i, simplify = FALSE)
  for (combo in combinations) {
    # Always include paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b, "_T1IFN_expression")
    covariates <- c(combo, constant_covariate)
    # Run Cox analysis
    cox_ph_result <- get_cox(cox_input_df, "OS.from.enrollment..months.", "Dead", covariates)
    if (!is.null(cox_ph_result$cox)) {
      temp <- summary(cox_ph_result$cox)
      hr_list <- c(hr_list, temp$coefficients[constant_covariate, "exp(coef)"])
      p_val_list <- c(p_val_list, temp$coefficients[constant_covariate, "Pr(>|z|)"])
      hr_upper_list <- c(hr_upper_list, temp$conf.int[constant_covariate, "upper .95"])
      hr_lower_list <- c(hr_lower_list, temp$conf.int[constant_covariate, "lower .95"])
      concordance_list <- c(concordance_list, temp$concordance[["C"]])
      covariates_combination_list <- c(covariates_combination_list, paste(covariates, collapse = ", "))
    }
  }
}



df <- data.frame(covariate_combination = covariates_combination_list,
                 p_value = p_val_list,
                 hazard_ratio = hr_list,
                 hazard_ratio_upper_CI = hr_upper_list,
                 hazard_ratio_lower_CI = hr_lower_list,
                 concordance = concordance_list)
