library(ggplot2)
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)




source("clonal_replacement_calculation.R")
source("summary_plot.R")

setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/"
)



get_clonal_replacement_with_survival_df = function(clonal_replacement_df,
                                              patient_survival_data,
                                              top)
{
  comparison_df = data.frame(Patient = rownames(clonal_replacement_df),
                             comparison_value = clonal_replacement_df[, paste0("top_", top)])
  select_columns =  c(
    "Patient",
    "PFS.from.enrollment..months." ,
    "OS.from.enrollment..months." ,
    "Dead",
    "Age",
    "Sex",
    "MGMT.methylation",
    "IDH.1.mutation",
    "Extension.of.resection..Biopsy.only..Partial..GTR."
    
  )
  if (is.null(comparison_df))
  {
    return(NULL)
  } else
  {
    comparison_df = merge(comparison_df,
                          patient_survival_data[, select_columns],
                          by = "Patient",
                          sort = FALSE)
    colnames(comparison_df)[c(3, 4)] = c("PFS", "OS")
    comparison_df$Dead = sapply(comparison_df$Dead, to_bool)
    comparison_df = comparison_df %>% filter(!is.na(comparison_value), !is.infinite(comparison_value))
    colnames(comparison_df)[ncol(comparison_df)]="resection"
    comparison_df$Max_resection = comparison_df$resection != "Biopsy"
    
    return(comparison_df)
  }
}



make_clonal_replacement_plot = function(comparison_df,
                                   compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                   plot_folder = "outs",
                                   variable_name = "clonal replacement",
                                   filename, clonal_replacement_method = "division")
{
  #Clonal_replacement_method: "substration" or division
  #alternative	: used in t test function. a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
  test_values = comparison_df$comparison_value
  test_values[test_values > 5] = 5 # Not proper test as test_values are always more than 1. Need to define variance somehow.
  print(paste("clonal replacement method",  clonal_replacement_method))
  if (clonal_replacement_method == "substraction")
  {
    mu  = 0
  }else
    if (clonal_replacement_method == "division")
    {
      mu  = 1
    }
  
  
  p_val = t.test(x = test_values,
                 mu = mu,
                 alternative = "g")[["p.value"]]
  p_label = paste("p =", round(p_val, 3))
  if (p_val <= 0.05)
  {
    col = "red"
    p_label = paste(p_label, "*")
    
  } else
  {
    col = "black"
  }
  
  
  plot_label = paste(variable_name ,
                     "\n",
                     compared_pair[[1]] ,
                     "vs",
                     compared_pair[[2]])
  x_label = paste(variable_name)
  p = ggplot(comparison_df, aes(x = x_label , y = comparison_value)) +
    geom_violin(trim = FALSE, width = 0.2) +
    geom_boxplot(width = 0.1) +
    geom_point() +
    annotate(
      "text",
      label = p_label,
      x = 1,
      y = min(max(comparison_df$comparison_value) * 1.5, 2.8),
      size = 20  ,
      col = col
    ) +
    scale_y_continuous(name = variable_name, limits = c(0,3)) +
    ggtitle(plot_label)  +
    geom_hline(
      yintercept = mu,
      linetype = "dashed",
      color = "red",
      linewidth = 2
    )
  
  p = plot_formating(p)  + theme(legend.position = "none")
  
  save_plot(p, filename = filename)
  return(p)
  
}



get_summarized_plots_clonal_replacement = function(comparison_df,
                                              patients,
                                              metadata,
                                              patient_survival_data,
                                              variable_name = "clonal replacement",
                                              compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                              plot_folder = "outs",
                                              clinical_formula =  'Age + Sex + MGMT.methylation + IDH.1.mutation')
{
  filename = paste0(plot_folder,
                    "/",
                    variable_name,
                    compared_pair[[1]],
                    "_vs_",
                    compared_pair[[2]],
                    ".pdf")
  clonal_replacement_plot = tryCatch({
    make_clonal_replacement_plot(
      comparison_df,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      variable_name = variable_name,
      filename = filename
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
      variable_name = variable_name,
      clinical_formula = clinical_formula
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
      variable_name = variable_name,
      clinical_formula = clinical_formula
    )
  },
  error = function(x)
  {
    print(x)
    return(NULL)
  })
  
  
  
  
  #combined_plots_list = c(list(clonal_replacement_plot), PFS_plots, OS_plots)
  
  combined_plots_list = c(list(clonal_replacement_plot),OS_plots)
  combined_plots = plot_grid(
    plotlist = combined_plots_list,
    nrow = 1,
    ncol = length(combined_plots_list),
    rel_widths = rep(1, length(combined_plots_list)),
    align = "h"
  )
  
  
  #
  # filename = paste0(
  #   plot_folder,
  #   "/",
  #   variable_name,"_",
  #   compared_pair[[1]],
  #   "_vs_",
  #   compared_pair[[2]],
  #   "_plots.pdf"
  # )
  # filename = gsub(" ", "_", filename)
  #
  # plot_width = 20 * 7
  # plot_height = 20
  # cairo_pdf(file =  filename ,
  #           width = plot_width,
  #           height = plot_height)
  # print(combined_plots)
  # dev.off()
  return(combined_plots)
}








make_top_plots = function(tops  = c(10, 20, 50, 100),
                          compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                          trackClonotype_df,
                          patient_survival_data,
                          metadata, plot_folder = "outs",
                          method = "division",
                          clinical_formula ='Age + Sex + MGMT.methylation + IDH.1.mutation')
{
  clonal_replacement_list =  get_top_sample_clones_change_clonal_timepoints(
    compared_pair = compared_pair,
    metadata = metadata,
    tops  = tops,
    trackClonotype_df = trackClonotype_df,
    method =  method
  )
  clonal_replacement_df = clonal_replacement_list$clonal_replacement_df
  top_plot_list = list()
  
  for (top in tops)
    
  {
    print(top)
    variable_name = paste("top", top, "clonal replacement")
    comparison_df = get_clonal_replacement_with_survival_df(
      clonal_replacement_df = clonal_replacement_df,
      patient_survival_data = patient_survival_data,
      top = top
    )

    
    p = get_summarized_plots_clonal_replacement(
      comparison_df = comparison_df,
      patients = patients,
      metadata = metadata,
      patient_survival_data = patient_survival_data,
      variable_name = variable_name,
      compared_pair = compared_pair,
      plot_folder = plot_folder,
      clinical_formula = clinical_formula
    )
    
    title <- ggdraw() +
      draw_label(
        paste("Top", top, " clonal replacement", compared_pair[[1]], compared_pair[[2]]),
        fontface = 'bold',
        x = 0,
        hjust = 0,
        size = 120
      ) +
      theme(# add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7))
    combined_plots = plot_grid(title, p,
                               ncol = 1,
                               # rel_heights values control vertical title margins
                               rel_heights = c(0.1, 1))
    
    combined_plots =     combined_plots + theme(# add margin on the left of the drawing canvas
      plot.margin = margin(0, 0, 0, 20))
    
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
    
    plot_width = 22 * length(p)
    plot_height = 20
    cairo_pdf(file =  filename ,
              width = plot_width,
              height = plot_height)
    print(combined_plots)
    dev.off()
    top_plot_list[[paste0("top", top)]] = combined_plots
  }
  
  top_plots = plot_grid(plotlist = top_plot_list, nrow = 1)
  
  filename = paste0(
    plot_folder,
    "/",
    "all_tops_clonal_replacement_",
    compared_pair[[1]],
    "_vs_",
    compared_pair[[2]],
    "_plots.pdf"
  )
  filename = gsub(" ", "_", filename)
  
  plot_width = 20 * 7 * length(tops)
  plot_height = 20
  cairo_pdf(file =  filename ,
            width = plot_width,
            height = plot_height)
  print(top_plots)
  dev.off()
  
  return(top_plots)
  
  
}


run_example = function()
  
  
{
trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt"

metadata = get_metadata_sorted()
trackClonotype_df <-
  read.delim(trackClonotype_df_file, check.names = FALSE)
patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
patient_survival_data <- read.delim(patient_survival_data_file)
patient_survival_data$Patient = paste0("p", patient_survival_data$Subject.ID)

make_top_plots(
  tops  = c(10, 20, 50, 100),
  compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
  trackClonotype_df = trackClonotype_df,
  patient_survival_data = patient_survival_data,
  metadata = metadata
)

}

#run_example()




