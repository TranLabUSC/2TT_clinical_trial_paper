# TCR clonal Replacement - Biopsy vs Resection Comparison
# Timepoints N_1, N_2, Y_2

library(dplyr)
library(tidyr)
library(ggplot2)

for (n in c(10, 20, 50, 100)){
  top_n = paste0("top_", n)
  
  survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
  survival_df <- read.csv(survival_data, sep = "\t")
  
  timepoint_a <- "Pre-TTF"
  timepoint_b <- "1"
  
  # read cloanl replacement ratio file
  clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
  clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
  clonal_replacement_ratio <- clonal_replacement_ratio[, top_n, drop = FALSE]
  names(clonal_replacement_ratio)[names(clonal_replacement_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
  rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
  
  survival_df <- merge(survival_df, clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
  
  
  timepoint_a <- "Pre-TTF"
  timepoint_b <- "2"
  
  # read cloanl replacement ratio file
  clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
  clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
  clonal_replacement_ratio <- clonal_replacement_ratio[, top_n, drop = FALSE]
  names(clonal_replacement_ratio)[names(clonal_replacement_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
  rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
  
  survival_df <- merge(survival_df, clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
  
  timepoint_a <- "Post-TTF"
  timepoint_b <- "2"
  
  # read cloanl replacement ratio file
  clonal_replacement_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD4_clonal_analysis/clonal_replacement_T_cells_division/"
  clonal_replacement_ratio <- read.table(paste0(clonal_replacement_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
  clonal_replacement_ratio <- clonal_replacement_ratio[, top_n, drop = FALSE]
  names(clonal_replacement_ratio)[names(clonal_replacement_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
  rownames(clonal_replacement_ratio) <- as.numeric(sub("^p", "", rownames(clonal_replacement_ratio)))
  
  survival_df <- merge(survival_df, clonal_replacement_ratio, by.x = "Subject.ID", by.y = "row.names")
  
  names(survival_df)[names(survival_df) == "Extension.of.resection..Biopsy.only..Partial..GTR."] <- "resection_type"
  values_to_replace <- c("GTR", "Partial")
  new_value <- "Maximal Resection"  
  survival_df[, "resection_type"][survival_df[, "resection_type"] %in% values_to_replace] <- new_value
  values_to_replace <- c("Biopsy")
  new_value <- "Biopsy Only"
  survival_df[, "resection_type"][survival_df[, "resection_type"] %in% values_to_replace] <- new_value
  
  
  # Convert the DataFrame to long format
  df_long <- pivot_longer(survival_df, cols = starts_with("Blood"), names_to = "count_type", values_to = "count_value")
  df_long$count_type <- factor(df_long$count_type, levels = c("Blood_Pre-TTF_vs_Blood_1", "Blood_Pre-TTF_vs_Blood_2", "Blood_Post-TTF_vs_Blood_2"))
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
    labs(title = paste0(top_n, " Clonal Replacement by Resection Type"), x = "Timepoint Comparison", y = "Clonal Replacement Ratio") +
    theme_minimal() +
    ylim(0,3) + 
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
  file_name <- paste0("Clonal_Replacement_Ratio_", top_n, "_Resection_Type_Comparison.pdf")
  ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 14, height = 4, dpi = 300, limitsize = FALSE)
}

