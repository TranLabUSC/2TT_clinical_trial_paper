#Plot compare two time points: preTTF vs all other time point
#post TTF vs all time point.
library(dplyr)
library(data.table)
library(ggplot2)
library(survival)
library(ggfortify)
library("survminer")
library("Rcpp")
library(cowplot)
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("get_metadata.R")




comparison_function = function(x, y)
{
  return(y / x)
}



get_comparison_df = function(value_to_compare_df,
                             compared_pair,
                             comparison_function,
                             patient_survival_data,
                             patients,
                             metadata)
{
  comparison_df = NULL
  for (patient in patients)
  {
    sample1 = metadata %>% filter(Patient == patient,
                                  integrated_sample_type == compared_pair[[1]])
    sample1 = sample1$sample_id
    
    sample2 = metadata %>% filter(Patient == patient,
                                  integrated_sample_type == compared_pair[[2]]) %>% select(sample_id)
    sample2 = sample2$sample_id
    if (length(sample1) > 0  & length(sample2) > 0)
    {
      value1 = value_to_compare_df[sample1, 1]
      value2 = value_to_compare_df[sample2, 1]
      pair_df = data.frame(
        Patient = patient,
        t1 = value1,
        t2 = value2 ,
        comparison_value = comparison_function(value1, value2)
      )
      if (is.null(comparison_df))
      {
        comparison_df = pair_df
      } else
      {
        comparison_df = rbind(comparison_df, pair_df)
      }
    }
    
  }
  select_columns =  c(
    "Patient",
    "PFS.from.enrollment..months." ,
    "OS.from.enrollment..months." ,
    "Dead",
    "Age",
    "Sex",
    "MGMT.methylation",
    "IDH.1.mutation"
    
  )
  if (is.null(comparison_df))
  {
    return(NULL)
  }else
  {  
  comparison_df = merge(comparison_df,
                        patient_survival_data[, select_columns],
                        by = "Patient",
                        sort = FALSE)
  colnames(comparison_df)[c(5, 6)] = c("PFS", "OS")
  comparison_df$Dead = sapply(comparison_df$Dead, to_bool)
  
  return(comparison_df)
  }
}

save_plot = function(plot, filename)
{
  filename = gsub(" ", "_", filename)
  plot_width = 15
  plot_height = 15
  cairo_pdf(file =  filename ,
            width = plot_width,
            height = plot_height)
  print(plot)
  dev.off()
  #return(p1)
}
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
      #axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 40)
    ) +
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 40
    ))
  
  return(p)
}


make_pairwise_plot = function(comparison_df,
                              compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                              plot_folder = "outs",
                              alterative = "g",
                              variable_name = "Shannon diversity")
  #Point pairwise comparison
{
  p_val = t.test(
    x = comparison_df$t1,
    y = comparison_df$t2,
    paired = TRUE,
    alternative = alterative
  )[["p.value"]]
  p_label = paste("p=", round(p_val, 3))
  if (p_val <= 0.05)
  {
    col = "red"
    p_label = paste(p_label, "*")
    
  } else
  {
    col = "black"
  }
  
  plot_label = paste(variable_name , "\n",  compared_pair[[1]] , "vs", compared_pair[[2]])
  timepoint_values_df = melt(setDT(comparison_df[, c("Patient", "t1", "t2")]),
                             id.vars = c("Patient"),
                             variable.name = "timepoint")
  timepoint_values_df  = as.data.frame(timepoint_values_df)
  max_val = max(timepoint_values_df$value)
  pairwise_plot = ggplot(timepoint_values_df, aes(x = timepoint, y = value)) +
    geom_violin(trim = FALSE, width = 0.2) +
    geom_boxplot(width = 0.1) +
    geom_point() +
    geom_line(aes(group = Patient, col = Patient))  +
    annotate(
      "text",
      label = p_label,
      x = 1.5,
      y = max_val * 1.1,
      size = 20  ,
      col = col
    ) +
    scale_x_discrete(labels = compared_pair) +
    scale_y_continuous(name = variable_name) +
    ggtitle(plot_label)
  
  
  # Formatting
  p = plot_formating(pairwise_plot)
  filename = paste0(
    plot_folder,
    "/",
    variable_name,
    "_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    ".pdf"
  )
  save_plot(p, filename = filename)
  return(p)
  
}

make_comparison_value_plot = function(comparison_df,
                                      compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                      plot_folder = "outs",
                                      alternative = "l",
                                      variable_name = "Shannon diversity")
{
  #alternative	: used in t test function. a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
  p_val = t.test(x = comparison_df$comparison_value,
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
                     variable_name ,
                     "\n",
                     compared_pair[[1]] ,
                     "vs",
                     compared_pair[[2]])
  x_label = paste("relative", variable_name)
  p = ggplot(comparison_df, aes(x = x_label , y = comparison_value)) +
    geom_violin(trim = FALSE, width = 0.2) +
    geom_boxplot(width = 0.1) +
    geom_point() +
    annotate(
      "text",
      label = p_label,
      x = 1,
      y = max(comparison_df$comparison_value) * 1.5,
      size = 20  ,
      col = col
    ) +
    scale_y_continuous(name = variable_name) +
    ggtitle(plot_label)  +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "red",
      linewidth = 2
    )
  
  p = plot_formating(p)  + theme(legend.position = "none")
  filename = paste0(
    plot_folder,
    "/",
    "relative_",
    variable_name,
    "_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    ".pdf"
  )
  save_plot(p, filename = filename)
  return(p)
  
}



#p1 = make_pairwise_plot(comparison_df = comparison_df)
#p2 = make_comparison_value_plot(comparison_df = comparison_df)
#print(p2)

#pfs scatter plot

make_survival_vs_comparison_value_plot = function(cox,
                                                  comparison_df,
                                                  survival_measure = "PFS",
                                                  compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                                  plot_folder = "outs",
                                                  variable_name = "Shannon diversity")
{
  coefs = summary(cox)$coefficients
  HR = coefs["comparison_value", "exp(coef)"]
  p_val =  coefs["comparison_value", "Pr(>|z|)"]
  p_label = paste(" CoxPH HR=", round(HR, 2), ", p=", round(p_val, 3))
  if (p_val <= 0.05)
  {
    col = "red"
    p_label = paste(p_label, "*")
    
  } else
  {
    col = "black"
  }
  
  plot_label = paste(
    "relative ",
    variable_name,
    survival_measure,
    " \n",
    compared_pair[[1]],
    "vs" ,
    compared_pair[[2]]
  )
  
  p = ggplot(data = comparison_df, aes(x = comparison_value, y = .data[[survival_measure]])) +   geom_point(aes(col =
                                                                                                                  Dead), size = 10) +
    annotate(
      "text",
      label = p_label,
      x = min(max(comparison_df$comparison_value) * 0.5, 2.5),
      y = max(comparison_df[, survival_measure]) * 0.95,
      size = 20  ,
      col = col
    ) + scale_x_continuous(limits = c(0, min(5, max(comparison_df$comparison_value))  )) +
    ggtitle(label = plot_label) + xlab(paste("relative", variable_name))
  
  
  p = plot_formating(p)
  
  
  
  filename = paste0(
    plot_folder,
    "/",
    survival_measure,
    "_vs_relative_",
    variable_name,
    "_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    ".pdf"
  )
  save_plot(p, filename = filename)
  return(p)
}


make_cox_fit = function(cox,
                        filename,
                        survival_measure,
                        compared_pair)
{
  cox_fit <- survfit(cox)
  p_cox_fit = autoplot(cox_fit)
  concordance = summary(cox)$concordance [["C"]]
  plot_label = paste(
    survival_measure,
    " Cox fit model\n",
    compared_pair[[1]],
    "vs" ,
    compared_pair[[2]],
    "\n Concordance=",
    round(concordance, 3)
  )
  p_cox_fit = p_cox_fit + ggtitle(plot_label) + xlab(survival_measure)
  p_cox_fit = plot_formating(p_cox_fit)
  
  
  save_plot(p_cox_fit, filename = filename)
  
  return(p_cox_fit)
}

make_comparison_value_low_high_fit_plot = function(comparison_df,
                                                   survival_measure = "PFS",
                                                   compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                                   plot_folder = "outs",
                                                   variable_name = "Shannon diversity")
{
  comparison_value_median = median(comparison_df$comparison_value)
  comparison_df$comparison_value_group = sapply(comparison_df$comparison_value, function(x)
  {
    if (x <= comparison_value_median)
    {
      return("Low")
    }
    else{
      return("High")
    }
  })
  
  survival_object = as.formula(
    paste0(
      "Surv(comparison_df[, '",
      survival_measure,
      "'], comparison_df$Dead)  ~   comparison_value_group"
    )
  )
  
  km_comparison_value_group_fit <-
    survfit(survival_object, data = comparison_df)
  surv_diff <- survdiff(survival_object, data = comparison_df)
  
  p_val = round(pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = FALSE), 3)
  #p_val = round(surv_diff$pvalue, 3)
  p_label = paste0("logrank p-value=", p_val)
  if (p_val <= 0.05)
  {
    col = "red"
    p_label = paste(p_label, "*")
    
  } else
  {
    col = "black"
  }
  
  p = autoplot(km_comparison_value_group_fit)
  p = p +   annotate(
    "text",
    label = p_label,
    y = 0.95,
    x = max(comparison_df[, survival_measure]) * 0.5,
    size = 20  ,
    col = col
  )
  
  plot_label = paste(
    survival_measure,
    " KM  plot\n stratified by relative ",
    variable_name,
    " \n",
    compared_pair[[1]],
    "vs" ,
    compared_pair[[2]]
  )
  p = p + ggtitle(plot_label) + xlab(survival_measure)
  p = plot_formating(p)
  
  filename = paste0(
    plot_folder,
    "/",
    survival_measure,
    "_KM_plot",
    
    "_by_relative_",
    variable_name,
    "_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    ".pdf"
  )
  save_plot(p, filename = filename)
  
  return(p)
}

make_survival_plots = function(comparison_df ,
                               survival_measure = "PFS",
                               compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                               plot_folder = "outs",
                               variable_name = "Shannon diversity")
{
  survival_object = as.formula(
    paste0(
      "Surv(comparison_df[, '",
      survival_measure,
      "'], comparison_df$Dead)  ~  Age + Sex + MGMT.methylation + IDH.1.mutation + comparison_value"
    )
  )
  cox  = coxph(survival_object, data = comparison_df, na.action = na.omit)
  
  
  survival_dot_plot = tryCatch({
     make_survival_vs_comparison_value_plot(
      cox = cox,
      comparison_df = comparison_df ,
      survival_measure = survival_measure,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      variable_name = variable_name
    )
    
    
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  
  cox_fit_filename = paste0(
    plot_folder,
    "/",
    "coxfit_",
    survival_measure,
    "_vs_relative_",
    variable_name,
    "_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    ".pdf"
  )
  
  cox_fit_plot = tryCatch({
    make_cox_fit(
      cox = cox,
      filename = cox_fit_filename,
      survival_measure = survival_measure,
      compared_pair = compared_pair
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  
  km_plot = tryCatch({
    make_comparison_value_low_high_fit_plot(
      comparison_df = comparison_df ,
      survival_measure = survival_measure,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      variable_name = variable_name
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  plot_list = list(survival_dot_plot, cox_fit_plot, km_plot)
  survival_plots = plot_grid(
    plotlist = plot_list,
    nrow = 1,
    ncol = 3,
    align = "h"
  )
  
  filename = paste0(plot_folder, "/", survival_measure, "_plots.pdf")
  plot_width = 20 * 3
  plot_height = 20
  cairo_pdf(file =  filename ,
            width = plot_width,
            height = plot_height)
  print(survival_plots)
  dev.off()
  return(plot_list)
}


get_summarized_plots = function(value_to_compare_df,
                                patients,
                                metadata,
                                patient_survival_data,
                                variable_name = "Shannon diversity",
                                compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                plot_folder = "outs",
                                pairwise_plot_alternative = "g",
                                relative_plot_alternative = "l")
{
  comparison_df = get_comparison_df(
    value_to_compare_df = value_to_compare_df,
    compared_pair = compared_pair,
    comparison_function = comparison_function,
    patient_survival_data = patient_survival_data,
    patients = patients,
    metadata = metadata
  )
  
  pairwise_plot = tryCatch({
    make_pairwise_plot(
      comparison_df,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      alterative = pairwise_plot_alternative,
      variable_name = variable_name
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  relative_plot = tryCatch({
    make_comparison_value_plot(
      comparison_df,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      alternative = relative_plot_alternative,
      variable_name = variable_name
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  
  
  PFS_plots = tryCatch({
    make_survival_plots(
      comparison_df,
      survival_measure = "PFS",
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      variable_name = variable_name
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  
  OS_plots = tryCatch({
    make_survival_plots(
      comparison_df,
      survival_measure = "OS" ,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      variable_name = variable_name
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  
  combined_plots_list = c(list(pairwise_plot), list(relative_plot), PFS_plots, OS_plots)
  combined_plots = plot_grid(
    plotlist = combined_plots_list,
    nrow = 1,
    ncol = 8,
    rel_widths = rep(1,8),
    align = "h"
  )
  
  filename = paste0(
    plot_folder,
    "/",
    variable_name,
    "_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    "_plots.pdf"
  )
  filename = gsub(" ", "_", filename)
  
  plot_width = 20 * 8
  plot_height = 20
  cairo_pdf(file =  filename ,
            width = plot_width,
            height = plot_height)
  print(combined_plots)
  dev.off()
  return(combined_plots)
}
#return(p1)


run_example = function()
{
  
  setwd(
    "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
  )
  
  shannon_clonal_diversity <-
    read.csv(
      "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/shannon_clonal_diversity.txt",
      sep = ""
    )
  metadata = get_metadata_sorted()
  patients = unique(metadata$Patient)
  patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
  patient_survival_data <- read.delim(patient_survival_data_file)
  patient_survival_data$Patient = paste0("p", patient_survival_data$Subject.ID)
  
  get_summarized_plots(
    value_to_compare_df = shannon_clonal_diversity,
    patients = patients,
    metadata = metadata,
    patient_survival_data = patient_survival_data,
    variable_name = "Shannon diversity",
    compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
    plot_folder = "outs",
    pairwise_plot_alternative = "g",
    relative_plot_alternative = "l"
  )
}

#run_example()
