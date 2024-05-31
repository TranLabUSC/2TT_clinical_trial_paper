library(ggplot2)
clonotyping_dir = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
setwd(paste0(clonotyping_dir, "/code"))
source("clonal_replacement_plot.R")


setwd(clonotyping_dir)








make_all_comparisons_plots = function(compared_pair_list,
                                      metadata,
                                      trackClonotype_df,
                                      patient_survival_data,
                                      plot_folder = "outs",
                                      tops = c(10, 20, 50, 100),
                                      method = "division",
                                      clinical_formula = 'Age + Sex + MGMT.methylation + IDH.1.mutation')
{
  dir.create(plot_folder)
  
  plot_list = list()
  for (compared_pair in compared_pair_list)
  {
    print(
      paste(
        "Making clonal replacement comparison",
        compared_pair[[1]],
        "-",
        compared_pair[[2]]
      )
    )
    print(paste("clinical formula :", clinical_formula))
    p = make_top_plots(
      tops  = tops,
      compared_pair = compared_pair,
      trackClonotype_df = trackClonotype_df,
      patient_survival_data = patient_survival_data,
      metadata = metadata,
      plot_folder = plot_folder,
      method = method,
      clinical_formula = clinical_formula
    )
    
    
    plot_list[[paste0(compared_pair[[1]], "_", compared_pair[[2]])]] = p
    
  }
  
  all_comparisons_plot = plot_grid(plotlist = plot_list,
                                   ncol = 1,
                                   nrow = length(plot_list))
  
  filename = paste0(plot_folder, "/clonal_replacement_timepoints_plots.pdf")
  filename = gsub(" ", "_", filename)
  
  plot_width = 20 * 8 * length(tops)
  plot_height = 25 * length(plot_list)
  cairo_pdf(file =  filename ,
            width = plot_width,
            height = plot_height)
  print(all_comparisons_plot)
  dev.off()
  return(all_comparisons_plot)
  
}



get_timepoint_clonal_replacement_plot = function(compared_pair_list,
                                                 trackClonotype_df_file = paste0(clonotyping_dir, "/trackClonotype_df.proportion.txt"),
                                                 patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv",
                                                 GBM_only = TRUE,
                                                 plot_folder = "outs_testing_top10_GBM_only",
                                                 tops = c(10),
                                                 method = "division",
                                                 clinical_formula_intent_to_treat = 'Age + Sex + MGMT.methylation + IDH.1.mutation',
                                                 clinical_formula_GBM_only =  'Age + Sex + MGMT.methylation')
  
{
  metadata = get_metadata_sorted()
  trackClonotype_df <-
    read.delim(trackClonotype_df_file, check.names = FALSE)
  
  patient_survival_data <- read.delim(patient_survival_data_file)
  patient_survival_data$Patient = paste0("p", patient_survival_data$Subject.ID)
  
  if (GBM_only)
  {
    patient_survival_data = patient_survival_data[patient_survival_data$IDH.1.mutation ==
                                                    "negative", ]
    clinical_formula = clinical_formula_GBM_only
    
  } else
  {
    clinical_formula =  clinical_formula_intent_to_treat
  }
  
  
  
  make_all_comparisons_plots(
    compared_pair_list = compared_pair_list,
    metadata = metadata,
    trackClonotype_df = trackClonotype_df,
    patient_survival_data = patient_survival_data,
    plot_folder = plot_folder,
    tops = tops,
    method = method,
    clinical_formula = clinical_formula
  )
  
  
}


get_timepoint_clonal_replacement_plot_both_patient_group = function(
  compared_pair_list ,
  trackClonotype_df_file = paste0(clonotyping_dir, "/trackClonotype_df.proportion.txt"),
  patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv",
  plot_folder = "outs_GBM_only_ITT",
  tops = c(10),
  method = "division",
  clinical_formula_intent_to_treat = 'Age + Sex + MGMT.methylation + IDH.1.mutation',
  clinical_formula_GBM_only =  'Age + Sex + MGMT.methylation'
)

  
{
 dir.create(plot_folder)
   for (GBM_only in c(TRUE, FALSE))
  {
     if (GBM_only == TRUE)
     {
       new_plot_folder = paste0(plot_folder, "/GBM_Only")
     }else
     {
       new_plot_folder = paste0(plot_folder, "/ITT")
     }
     dir.create(new_plot_folder)
  get_timepoint_clonal_replacement_plot(
    compared_pair_list ,
    trackClonotype_df_file ,
    patient_survival_data_file,
    GBM_only = GBM_only,
    plot_folder = new_plot_folder,
    tops = tops,
    method = method,
    clinical_formula_intent_to_treat = clinical_formula_intent_to_treat,
    clinical_formula_GBM_only = clinical_formula_GBM_only
  )

   }

}



# get_timepoint_clonal_replacement_plot(
#   compared_pair_list = compared_pair_list,
#   trackClonotype_df_file = paste0(clonotyping_dir, "/trackClonotype_df.proportion.txt"),
#   patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv",
#   GBM_only = TRUE,
#   plot_folder = "outs_both_GBM_and_ITT",
#   tops = c(10),
#   method = "division"
# )


compared_pair_list = list(
  c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
  c("Blood_Pre-TTF_sc", "Blood_1_sc"),
  c("Blood_Pre-TTF_sc", "Blood_2_sc"),
  #c("Blood_Pre-TTF_sc", "Blood_3_sc" ),
  #c("Blood_Pre-TTF_sc", "Blood_4_sc" ),
  #c("Blood_Pre-TTF_sc", "Blood_5_sc" ),
  #c("Blood_Pre-TTF_sc", "Blood_R1_sc" ),
  #c("Blood_Pre-TTF_sc", "Blood_R2_sc" ),
  
  c("Blood_Post-TTF_sc", "Blood_1_sc"),
  c("Blood_Post-TTF_sc", "Blood_2_sc")
  #c("Blood_Post-TTF_sc", "Blood_3_sc" ),
  #c("Blood_Post-TTF_sc", "Blood_4_sc" ),
  #c("Blood_Post-TTF_sc", "Blood_5_sc" ),
  #c("Blood_Post-TTF_sc", "Blood_R1_sc" ),
  #c("Blood_Post-TTF_sc", "Blood_R2_sc" ),
  
  #c("Blood_1_sc" , "Blood_2_sc" ),
  #c("Blood_2_sc" , "Blood_3_sc" ),
  #c("Blood_3_sc" , "Blood_4_sc" ),
  #c("Blood_4_sc" , "Blood_5_sc" ),
  #c("Blood_1_sc" , "Blood_R1_sc" ),
  #c("Blood_2_sc" , "Blood_R1_sc" ),
  #c("Blood_3_sc" , "Blood_R1_sc" ),
  #c("Blood_4_sc" , "Blood_R1_sc" ),
  #c("Blood_5_sc" , "Blood_R1_sc" )
)



# get_timepoint_clonal_replacement_plot_both_patient_group(
#     compared_pair_list = compared_pair_list,
#     trackClonotype_df_file = paste0(clonotyping_dir, "/trackClonotype_df.proportion.all_T_cells.Harshit.txt"),
#     patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv",
#     plot_folder = "outs_new_harshit_clonal_replacement",
#     tops = c(10),
#     method = "division",
#     clinical_formula_intent_to_treat = 'MGMT.methylation + IDH.1.mutation + Age + Sex ',
#     clinical_formula_GBM_only =  'MGMT.methylation  + Age + Sex'
# )
#   

get_timepoint_clonal_replacement_plot_both_patient_group(
  compared_pair_list = compared_pair_list,
  trackClonotype_df_file = paste0(clonotyping_dir, "/trackClonotype_df.proportion.all_T_cells.Harshit.txt"),
  patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv",
  plot_folder = "outs_new_harshit_clonal_replacement_no_color",
  tops = c(10),
  method = "division",
  clinical_formula_intent_to_treat = ' IDH.1.mutation + Age + Sex ',
  clinical_formula_GBM_only =  ' Age + Sex'
)


