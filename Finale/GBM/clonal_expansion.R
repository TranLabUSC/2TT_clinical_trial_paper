# TCR clonal Expansion - Biopsy vs Resection Comparison
# Timepoints N_Y, N_1, N_2, Y_1, Y_2

library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(broom)
library(ggfortify)
library("survminer")
library("Rcpp")
library(cowplot)
library(tidyr)
library(rlang)


survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")
survival_df <- survival_df[survival_df$IDH.1.mutation == "negative", ]
top_n = "clonal_expansion"

timepoint_a <- "Pre-TTF"
timepoint_b <- "Post-TTF"

# read cloanl expansion ratio file
clonal_expansion_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
clonal_expansion_ratio <- read.table(paste0(clonal_expansion_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
clonal_expansion_ratio <- clonal_expansion_ratio[, top_n, drop = FALSE]
names(clonal_expansion_ratio)[names(clonal_expansion_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
rownames(clonal_expansion_ratio) <- as.numeric(sub("^p", "", rownames(clonal_expansion_ratio)))

survival_df <- merge(survival_df, clonal_expansion_ratio, by.x = "Subject.ID", by.y = "row.names", all.x = TRUE)


timepoint_a <- "Pre-TTF"
timepoint_b <- "1"

# read cloanl expansion ratio file
clonal_expansion_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
clonal_expansion_ratio <- read.table(paste0(clonal_expansion_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
clonal_expansion_ratio <- clonal_expansion_ratio[, top_n, drop = FALSE]
names(clonal_expansion_ratio)[names(clonal_expansion_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
rownames(clonal_expansion_ratio) <- as.numeric(sub("^p", "", rownames(clonal_expansion_ratio)))

survival_df <- merge(survival_df, clonal_expansion_ratio, by.x = "Subject.ID", by.y = "row.names", all.x = TRUE)

timepoint_a <- "Pre-TTF"
timepoint_b <- "2"

# read cloanl expansion ratio file
clonal_expansion_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
clonal_expansion_ratio <- read.table(paste0(clonal_expansion_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
clonal_expansion_ratio <- clonal_expansion_ratio[, top_n, drop = FALSE]
names(clonal_expansion_ratio)[names(clonal_expansion_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
rownames(clonal_expansion_ratio) <- as.numeric(sub("^p", "", rownames(clonal_expansion_ratio)))

survival_df <- merge(survival_df, clonal_expansion_ratio, by.x = "Subject.ID", by.y = "row.names", all.x = TRUE)


timepoint_a <- "Post-TTF"
timepoint_b <- "1"

# read cloanl expansion ratio file
clonal_expansion_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
clonal_expansion_ratio <- read.table(paste0(clonal_expansion_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
clonal_expansion_ratio <- clonal_expansion_ratio[, top_n, drop = FALSE]
names(clonal_expansion_ratio)[names(clonal_expansion_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
rownames(clonal_expansion_ratio) <- as.numeric(sub("^p", "", rownames(clonal_expansion_ratio)))

survival_df <- merge(survival_df, clonal_expansion_ratio, by.x = "Subject.ID", by.y = "row.names", all.x = TRUE)


timepoint_a <- "Post-TTF"
timepoint_b <- "2"

# read cloanl expansion ratio file
clonal_expansion_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
clonal_expansion_ratio <- read.table(paste0(clonal_expansion_dir, "Blood_", timepoint_a, "_sc_vs_Blood_", timepoint_b, "_sc.txt"), sep = "\t")
clonal_expansion_ratio <- clonal_expansion_ratio[, top_n, drop = FALSE]
names(clonal_expansion_ratio)[names(clonal_expansion_ratio) == top_n] <- paste0("Blood_", timepoint_a, "_vs_Blood_", timepoint_b)
rownames(clonal_expansion_ratio) <- as.numeric(sub("^p", "", rownames(clonal_expansion_ratio)))

survival_df <- merge(survival_df, clonal_expansion_ratio, by.x = "Subject.ID", by.y = "row.names", all.x = TRUE)

names(survival_df)[names(survival_df) == "Extension.of.resection..Biopsy.only..Partial..GTR."] <- "resection_type"
values_to_replace <- c("GTR", "Partial")
new_value <- "Maximal Resection"  
survival_df[, "resection_type"][survival_df[, "resection_type"] %in% values_to_replace] <- new_value
values_to_replace <- c("Biopsy")
new_value <- "Biopsy Only"
survival_df[, "resection_type"][survival_df[, "resection_type"] %in% values_to_replace] <- new_value


# Convert the DataFrame to long format
df_long <- pivot_longer(survival_df, cols = starts_with("Blood"), names_to = "count_type", values_to = "count_value")
df_long$count_type <- factor(df_long$count_type, levels = c("Blood_Pre-TTF_vs_Blood_Post-TTF", "Blood_Pre-TTF_vs_Blood_1", "Blood_Pre-TTF_vs_Blood_2", "Blood_Post-TTF_vs_Blood_1", "Blood_Post-TTF_vs_Blood_2"))
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
  geom_jitter(position = position_dodge(width = 1.2)) +  # Correct position dodging for jGBMer
  facet_wrap(~count_type, scales = "free") +  # Create a separate plot for each count type
  labs(title = paste0(top_n, " Clonal Diversity by Resection Type"), x = "Timepoint Comparison", y = "Clonal Diversity Ratio") +
  theme_minimal() +
  # ylim(0,3) + 
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

plot_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/"
file_name <- paste0("simpson_Clonal_Diversity_Ratio_Resection_Type_Comparison.pdf")
ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 12, height = 8, dpi = 300, limitsize = FALSE)


plot_formating = function(p)
{
  p = p + theme_bw() +
    theme(
      #legend.position = "none" ,
      text = element_text(size = 40),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 40)
    ) +
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 40
    ))
  
  return(p)
}
  
make_comparison_value_plot = function(df,
                                      variable_name,
                                      diversity_index = "simpson",
                                      alternative = "l",
                                      plot_dir)
{
  #alternative	: used in t test function. a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
  p_val = t.test(x = df[, variable_name],
                 mu = 1,
                 alternative = alternative)[["p.value"]]
  p_label = paste("p =", round(p_val, 3))
  if (p_val <= 0.05)
  {
    col = "red"
    p_label = paste(p_label, "*")
    
  } else
  {
    col = "black"
  }
  plot_label = paste("Relative",
                     variable_name)
  x_label = paste("relative", variable_name)
  if (diversity_index == "shannon"){
    p = ggplot(df, aes(x = x_label , y = !!sym(variable_name))) +
      geom_violin(trim = FALSE, width = 0.2) +
      geom_boxplot(width = 0.1) +
      geom_point() +
      annotate(
        "text",
        label = p_label,
        x = 1,
        y = max(df[, variable_name], na.rm = TRUE) * 1.5,
        size = 10  ,
        col = col
      ) +
      scale_y_continuous(name = variable_name) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_label)  +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "red",
        linewidth = 2
      )
  } else {
    p = ggplot(df, aes(x = x_label , y = !!sym(variable_name))) +
      geom_violin(trim = FALSE, width = 0.2) +
      geom_boxplot(width = 0.1) +
      geom_point() +
      annotate(
        "text",
        label = p_label,
        x = 1,
        y = max(df[, variable_name], na.rm = TRUE) * 1.1,
        size = 10  ,
        col = col
      ) +
      scale_y_continuous(name = variable_name) +
      ylab("Clonal Diversity Ratio") +
      ggtitle(plot_label)  +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "red",
        linewidth = 2
      )
  }
  
  p = p + theme_bw() +
    theme(
      legend.position = "none" ,
      text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_blank()
    ) +
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 10
    ))
  print(p)
  file_name <- paste0(diversity_index, "_Clonal_Diversity_", variable_name, ".pdf")
  ggsave(paste0(plot_dir, file_name), p, device = "pdf", width = 4, height = 6, dpi = 300, limitsize = FALSE)
}

make_comparison_value_plot(survival_df, "Blood_Pre-TTF_vs_Blood_Post-TTF", diversity_index = "simpson", alternative = "l", plot_dir)
make_comparison_value_plot(survival_df, "Blood_Pre-TTF_vs_Blood_1", diversity_index = "simpson", alternative = "l", plot_dir)
make_comparison_value_plot(survival_df, "Blood_Pre-TTF_vs_Blood_2", diversity_index = "simpson", alternative = "l", plot_dir)
make_comparison_value_plot(survival_df, "Blood_Post-TTF_vs_Blood_1", diversity_index = "simpson", alternative = "l", plot_dir)
make_comparison_value_plot(survival_df, "Blood_Post-TTF_vs_Blood_2", diversity_index = "simpson", alternative = "l", plot_dir)


############################################################################################################
# survival analysis with clonal diversity
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
survival_df <- survival_df[survival_df$IDH.1.mutation == "negative",]

timepoint_a <- "Pre-TTF"
timepoint_b <- "Post-TTF"
top_n <- "clonal_expansion"

# read clonal expansion ratio file
in_dir <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/CD8_clonal_analysis/clonal_expansion_simpson_all_T_cells_division/"
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
# cox_input_df[, top_n] <- 1/cox_input_df[, top_n]

# Define other covariates excluding 'top_10' which is always included
other_covariates <- c("Age", "Sex", "MGMT.methylation", "Extension.of.resection..Biopsy.only..Partial..GTR.")

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


