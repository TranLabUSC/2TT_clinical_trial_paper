devtools::install_version("ggplot2", "3.3.3")

library(ggplot2)
library(cowplot)



# Goal: 1) what are changes under Keytrude treatment : Gene, pathway.
#2) What are signature of response (Good, bad)

# Combine plots: 
#1) Clinical data
#2) PFS, OS
#3) Tumor volume
#4) Diversity
#4) Clonal expansion
#5)clonotype tracking
#6) Gene expression heatmap and UMAP of : Interferon gamma pathway. T Cell activation genes list, T cell inhibition gene list
#7) Cell type proportion change:DC, B, cell, .. accross time
#8) Gene expression of top 100 clones expanded.
#9) DC activation markers
#10) STING pathway?

#) 9 Tumor mutation burden
# Gene epxpression in blood of T cell clones that in the primary tumor and  recurrent tumor vs the rest.

setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("get_metadata.R")
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)

patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"

metadata = get_metadata_sorted()
patient_survival_data <- read.delim(patient_survival_data_file)
patients = unique(metadata$Patient)
patient_info_plots = readRDS("patient_info_plots_all_patients.RDS")
tumor_size_plots = readRDS("tumor_size_plots_all_patients.RDS")

Simpson_diversity_plots = readRDS("Simpson_diversity_plots_all_patients.RDS")
Shannon_diversity_plots = readRDS("Shannon_diversity_plots_all_patients.RDS")
clonal_distribution_plots = readRDS("clonal_distribution_plots_all_patients.RDS")
clonal_plots_top10 = readRDS("top10_clonal_plot_all_patients.RDS")
clonal_plots_top20 = readRDS("top20_clonal_plot_all_patients.RDS")
clonal_plots_top100 = readRDS("top100_clonal_plot_all_patients.RDS")
cell_clone_numbers_plots = readRDS("cell_clone_numbers_plots.RDS")


# For each patient, build the plot columns with following plots:
#1) Clinical info table
#2) OS, PFS
#3)Tumor Volume
#3b) cell_clone_numbers_plots
#4) Clonal Diversity
#5) Clonal Expansion
#6) Clonal tracking  top of all samples
# 7) "Tumor_P_immunoseq",

#9) "Blood_Pre-TTF_sc",
#10)"Blood_Post-TTF_sc",
#11)"Blood_1_sc",
#12)"Blood_2_sc",
#13)"Blood_3_sc",
#14)"Blood_4_sc",
#15)"Blood_5_sc",
#16)"Blood_R_sc",
#)17)"Blood_R2_sc",
#18)"Tumor_R_immunoseq",
#) 8)"Tumor_Tumor_sc",
#19)"Tumor_R1_immunoseq",
#20)"Tumor_R2_immunoseq"


# Make plot list for each patient
make_patient_plot_list = function(patient="p16", plot_folder = ".")
{
  
  #1) Clinical info table
  clinical_table = patient_info_plots[[patient]][["clinical_table" ]]
  #2) OS, PFS
  survival = patient_info_plots[[patient]][["survival_plot" ]]
  #3)Tumor Volume
  tumor_size = tumor_size_plots[[patient]]
  cell_clone_numbers =  cell_clone_numbers_plots[[patient]][["combined_plot"]]
  #4) Clonal Diversity
  simpson_clonal_diversity = Simpson_diversity_plots[[patient]][["clonal_diversity_plot"]]
  #+ theme(axis.text.x = element_blank(), axis.title.x =element_blank(), axis.ticks.x = element_blank() )
  #5) Clonal Expansion
  simpson_clonal_expanstion =  Simpson_diversity_plots[[patient]][["clonal_expansion_plot"]] 
  
  shannon_clonal_diversity = Shannon_diversity_plots[[patient]][["clonal_diversity_plot"]]
  #+ theme(axis.text.x = element_blank(), axis.title.x =element_blank(), axis.ticks.x = element_blank() )
  #5) Clonal Expansion
  shannon_clonal_expanstion =  Shannon_diversity_plots[[patient]][["clonal_expansion_plot"]] 
  clonal_distribution =  clonal_distribution_plots[[patient]]
  
  get_clonal_plots = function(clonal_plots)
  {
  top_all_samples_tracking = clonal_plots[[patient]][["combined_top_sequences_plot"]] 
  clonal_plot_list = clonal_plots[[patient]][["clonal_plot_list"]]
  Tumor_P_immunoseq = clonal_plot_list[["Tumor_P_immunoseq"]] 
  
  Blood_PreTTF_sc = clonal_plot_list[["Blood_Pre-TTF_sc"]]  
  Blood_PostTTF_sc= clonal_plot_list[["Blood_Post-TTF_sc"]]  
  Blood_1_sc =  clonal_plot_list[["Blood_1_sc"]] 
  Blood_2_sc = clonal_plot_list[["Blood_2_sc"]] 
  Blood_3_sc = clonal_plot_list[["Blood_3_sc"]] 
  Blood_4_sc = clonal_plot_list[["Blood_4_sc"]] 
  Blood_5_sc = clonal_plot_list[["Blood_5_sc"]]
  Blood_R_sc = clonal_plot_list[["Blood_R_sc"]] 
  Blood_R2_sc = clonal_plot_list[["Blood_R2_sc"]] 
  Tumor_Tumor_sc = clonal_plot_list[["Tumor_Tumor_sc"]] 
  Tumor_R_immunoseq = clonal_plot_list[["Tumor_R_immunoseq"]]  
  Tumor_R1_immunoseq = clonal_plot_list[["Tumor_R1_immunoseq" ]] 
  Tumor_R2_immunoseq =  clonal_plot_list[["Tumor_R2_immunoseq" ]]  
  clonal_plot_list = list(
                          top_all_samples_tracking = top_all_samples_tracking ,
                          Tumor_P_immunoseq = Tumor_P_immunoseq ,
                          Blood_PreTTF_sc = Blood_PreTTF_sc ,
                          Blood_PostTTF_sc = Blood_PostTTF_sc,
                          Blood_1_sc =  Blood_1_sc ,
                          Blood_2_sc = Blood_2_sc,
                          Blood_3_sc = Blood_3_sc,
                          Blood_4_sc = Blood_4_sc,
                          Blood_5_sc = Blood_5_sc,
                          Blood_R_sc = Blood_R_sc,
                          Blood_R2_sc = Blood_R2_sc,
                          Tumor_Tumor_sc =  Tumor_Tumor_sc ,
                          Tumor_R_immunoseq = Tumor_R_immunoseq,
                          Tumor_R1_immunoseq = Tumor_R1_immunoseq,
                          Tumor_R2_immunoseq = Tumor_R2_immunoseq
  )
  
  clonal_plots = plot_grid(plotlist =   clonal_plot_list, ncol=1, nrow = length(clonal_plot_list)) 
  return(clonal_plots)
  }
  
  
  plot_list = list(clinical_table = clinical_table  ,
                   survival = survival ,
                   tumor_size = tumor_size ,
                   cell_clone_numbers = cell_clone_numbers,
                   simpson_clonal_diversity = simpson_clonal_diversity,
                   simpson_clonal_expanstion = simpson_clonal_expanstion,
                   shannon_clonal_diversity = shannon_clonal_diversity,
                   shannon_clonal_expanstion = shannon_clonal_expanstion,
                   clonal_distribution = clonal_distribution
                   )
      heading_plots =  plot_grid(plotlist = plot_list, ncol=1, rel_heights = c(1,0.4,1,4,1,1,1,1,13)) 
      
      top10_clonal_plots = get_clonal_plots(clonal_plots_top10)
      top20_clonal_plots = get_clonal_plots(clonal_plots_top20)
      top100_clonal_plots = get_clonal_plots(clonal_plots_top100)
      
      plot_list = list(heading_plots, top10_clonal_plots, top20_clonal_plots, top100_clonal_plots)
      all_plots = plot_grid(plotlist = plot_list, ncol=1, rel_heights =  c(23.4, 20,20,20)) 
      
      filename = paste0(plot_folder, "/", patient,  "_all_plots_combined.pdf")
      plot_width = 15
      plot_height = 700
      
      cairo_pdf(
        file =  filename ,
        width = plot_width,
        height = plot_height, pointsize = 24
      )
      print(all_plots)
      dev.off()
      print(paste("done combining patient plot for patient", patient))
      
      return(list(patient_plot_list = plot_list, patient_combined_plot = all_plots))

  
}
    



make_plots_for_all_patients = function(plot_folder ="outs")
{
  #unique(patient_survival_data$Best.ORR)
  #patient_survival_group = c("PD"  , "SD" ,  "PR", "PR/near CR" )
  #patient_survival_group = c("PR/near CR" )
  all_subplots = list()
  all_patient_plots = list()
  
  # 
  # 
  # for (grp in patient_survival_group)
  # {
  #   print(grp)
  #   grp_plot_list = list()
    #patient_survivals = patient_survival_data %>% filter(Best.ORR == grp) %>% arrange(OS.from.enrollment..months.)
    patient_survivals = patient_survival_data  %>% arrange(OS.from.enrollment..months.)
    patients = patient_survivals$Subject.ID
    patients = paste0("p", patients)
    for (patient in patients)
    {
      print(paste("making plot for patient", patient))
      plot_folder_patient = paste0(plot_folder, "/", patient)
      dir.create(plot_folder_patient)
       patient_plots = make_patient_plot_list(patient=patient, plot_folder =plot_folder_patient ) 
       patient_combined_plot  =  patient_plots[["patient_combined_plot" ]]
       patient_plot_list =  patient_plots[["patient_plot_list"]]
       #grp_plot_list[[patient]] = patient_combined_plot
       all_subplots[[patient]] = patient_plot_list
       all_patient_plots[[patient]] = patient_combined_plot
    }
    # Now print (group plot)
    # grp_combined_plot = plot_grid(plotlist=grp_plot_list, ncol=length(patients), nrow=1)
    # grp = gsub(" ","_", grp)
    # grp = gsub("/","_", grp)
    # 
    # filename = paste0(plot_folder, "/", grp,  ".pdf")
    # print(filename)
    # plot_width = 15 * length(grp_plot_list)
    # plot_height = 200
    # 
    # cairo_pdf(
    #   file =  filename ,
    #   width = plot_width,
    #   height = plot_height, pointsize = 24
    # )
  #   print(grp_combined_plot)
  #   dev.off()
  #   print(paste("done combining group plot for group", grp))
  #   
  # }
  
  # Now print (group plot)
  all_combined_plot = plot_grid(plotlist=all_patient_plots, nrow=1, ncol = length(all_patient_plots))
  filename = paste0(plot_folder, "/all_patients_plots.pdf")
  plot_width = 15 * length(all_patient_plots)
  plot_height = 700
  
  cairo_pdf(
    file =  filename ,
    width = plot_width,
    height = plot_height, pointsize = 24
  )
  print(all_combined_plot)
  dev.off()
  print("done combining all plots for all patients")
  return(list(all_patient_plots= all_patient_plots, combined_plot = all_combined_plot))
  
}


make_plots_for_all_patients_sorted_by_resection = function(plot_folder ="outs")
{
    
  dir.create(plot_folder)
  groups = c("Biopsy", "Partial", "GTR")
  group_plot_list  = list()
  for (group in groups)
  {
    all_subplots = list()
    all_patient_plots = list()
    patient_survivals = patient_survival_data %>% filter( Extension.of.resection..Biopsy.only..Partial..GTR. == group)  %>% arrange(OS.from.enrollment..months.)
    patients = patient_survivals$Subject.ID
    patients = paste0("p", patients)
    
    
    for (patient in patients)
    {
      print(paste("making plot for patient", patient))
      plot_folder_patient = paste0(plot_folder, "/", patient)
      dir.create(plot_folder_patient)
      patient_plots = make_patient_plot_list(patient=patient, plot_folder =plot_folder_patient ) 
      patient_combined_plot  =  patient_plots[["patient_combined_plot" ]]
      patient_plot_list =  patient_plots[["patient_plot_list"]]
      #grp_plot_list[[patient]] = patient_combined_plot
      all_subplots[[patient]] = patient_plot_list
      all_patient_plots[[patient]] = patient_combined_plot
    }
  
  
    all_combined_plot = plot_grid(plotlist=all_patient_plots, nrow=1, ncol = length(all_patient_plots))
    
    all_combined_plot <- draw_label(group, size = 20) + all_combined_plot
    group_plot_list[[group]] = all_combined_plot
    
    
    
    
    filename = paste0(plot_folder, "/", group, "_patients_plots.pdf")
    plot_width = 15 * length(all_patient_plots)
    plot_height = 700
    
    cairo_pdf(
      file =  filename ,
      width = plot_width,
      height = plot_height, pointsize = 24
    )
    print(all_combined_plot)
    dev.off()
    print(paste("done combining all plots for  patients group ", group))
  }
  
  
  all_groups_plot = plot_grid(plotlist=group_plot_list, nrow=1, ncol = 3)
  
  
  filename = paste0(plot_folder, "/all_groups_plots.pdf")
  plot_width = 15 * nrow(patient_survival_data)
  plot_height = 700
  
  cairo_pdf(
    file =  filename ,
    width = plot_width,
    height = plot_height, pointsize = 24
  )
  print(all_groups_plot)
  dev.off()
  print(paste("done combining all plots for  patients group ", group))
  
}



plot_folder ="outs_resection_group"
out = make_plots_for_all_patients_sorted_by_resection(plot_folder = plot_folder)

#saveRDS(out, file = "all_plots.RDS" )




