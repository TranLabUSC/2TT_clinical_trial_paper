setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("summary_plot.R")
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/"
)

sample_order =  c (
  "Tumor_P_immunoseq",
  "Blood_Pre-TTF_sc",
  "Blood_Post-TTF_sc",
  "Blood_1_sc",
  "Blood_2_sc",
  "Blood_3_sc",
  "Blood_4_sc",
  "Blood_5_sc",
  "Blood_R1_sc",
  "Blood_R2_sc",
  "Tumor_Tumor_sc",
  "Tumor_R1_immunoseq",
  "Tumor_R2_immunoseq"
)


metadata = get_metadata_sorted()
patients = unique(metadata$Patient)
patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
patient_survival_data <- read.delim(patient_survival_data_file)
patient_survival_data$Patient = paste0("p", patient_survival_data$Subject.ID)
compared_pair_list = list(
  c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
  c("Blood_Post-TTF_sc", "Blood_1_sc" ),
  c("Blood_Post-TTF_sc", "Blood_2_sc" ),
  c("Blood_Post-TTF_sc", "Blood_3_sc" ),
  c("Blood_Post-TTF_sc", "Blood_4_sc" ),
  c("Blood_Post-TTF_sc", "Blood_5_sc" ),
  c("Blood_Post-TTF_sc", "Blood_R1_sc" ),
  c("Blood_Post-TTF_sc", "Blood_R2_sc" ),
  c("Blood_1_sc" , "Blood_2_sc" ),
  c("Blood_2_sc" , "Blood_3_sc" ),
  c("Blood_3_sc" , "Blood_4_sc" ),
  c("Blood_4_sc" , "Blood_5_sc" ),
  c("Blood_1_sc" , "Blood_R1_sc" ),
  c("Blood_2_sc" , "Blood_R1_sc" ),
  c("Blood_3_sc" , "Blood_R1_sc" ),
  c("Blood_4_sc" , "Blood_R1_sc" ),
  c("Blood_5_sc" , "Blood_R1_sc" )
  
)
#compared_pair_list = list( c("Blood_5_sc" , "Blood_R1_sc" ))

make_all_comparisons_plots= function(value_to_compare_df = shannon_clonal_diversity,
                                     metadata,
                                     patients,
                                     variable_name = "Shannon diversity",
                                     plot_folder = "outs",
                                     pairwise_plot_alternative = "g",
                                     relative_plot_alternative = "l",
                                     compared_pair_list
                                     )
{
    plot_list = list()
  for (compared_pair in compared_pair_list)
  {
  print(paste("Making comparison", compared_pair[[1]], "-", compared_pair[[2]]))
  p =get_summarized_plots(value_to_compare_df = value_to_compare_df,
                       patients = patients,
                       metadata = metadata,
                       patient_survival_data = patient_survival_data,
                       variable_name = variable_name,
                       compared_pair = compared_pair,
                       plot_folder = plot_folder,
                       pairwise_plot_alternative = pairwise_plot_alternative,
                       relative_plot_alternative = relative_plot_alternative)
  
  plot_list[[paste0(compared_pair[[1]], "_", compared_pair[[2]])]] = p
  
  }
  
  all_comparisons_plot = plot_grid(plotlist = plot_list, ncol=1, nrow = length(plot_list))
  
  filename = paste0(plot_folder, "/", variable_name, "_timepoints_plots.pdf")
  filename = gsub(" ", "_", filename)
  
  plot_width = 20 * 8
  plot_height = 20 * length(plot_list)
  cairo_pdf(file =  filename ,
            width = plot_width,
            height = plot_height)
  print(all_comparisons_plot)
  dev.off()
  return(all_comparisons_plot)
  
}

get_shannon_diversity_timepoint_plots = function( diversity_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/shannon_clonal_diversity.txt",
                                                  plot_folder)
{
  shannon_clonal_diversity <-
    read.csv( diversity_file,sep = "")
      
make_all_comparisons_plots(
                        value_to_compare_df = shannon_clonal_diversity,
                            metadata = metadata,
                           patients = patients,
                           compared_pair_list = compared_pair_list,
                           variable_name = "Shannon diversity",
                        pairwise_plot_alternative = "g",
                        relative_plot_alternative = "l",
                        plot_folder = plot_folder
                        )
}


get_shannon_clonal_expansion_timepoint_plots = function( diversity_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/shannon_clonal_diversity.txt",
                                                         plot_folder)
{
  shannon_clonal_diversity <-
    read.csv( diversity_file,sep = "")
  shannon_clonal_expansion = 1/shannon_clonal_diversity
  colnames(shannon_clonal_expansion) = "clonal_expansion"
  
  make_all_comparisons_plots(
    value_to_compare_df = shannon_clonal_expansion,
    metadata = metadata,
    patients = patients,
    compared_pair_list = compared_pair_list,
    variable_name = "Shannon clonal expansion",
    pairwise_plot_alternative = "l",
    relative_plot_alternative = "g",
    plot_folder = plot_folder)
}




get_top_clonal_expansion_timepoint_plots = function(top = 10, top_proportions,
                                                    plot_folder)
{

  top_col = c(paste0("top_", top, "_proportion"))
  value_to_compare_df  = subset(top_proportions,select=c(top_col))
  make_all_comparisons_plots(
    value_to_compare_df = value_to_compare_df,
    metadata = metadata,
    patients = patients,
    compared_pair_list = compared_pair_list,
    variable_name = paste("top", top, "clonal expansion"),
    pairwise_plot_alternative = "l",
    relative_plot_alternative = "g",
    plot_folder = plot_folder)

}

get_tops_plots = function(top_proportion_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/top_proportions.txt",
                          tops = c(10,20,50,100),
                          plot_folder)
{
  top_proportions <- read.delim(top_proportion_file)
  for (top in tops)
  {
    print(top)
    get_top_clonal_expansion_timepoint_plots(top = top, top_proportions, plot_folder = plot_folder)
  }
  
}

#plot_folder = "outs"
#get_shannon_diversity_timepoint_plots(plot_folder = plot_folder)
#get_shannon_clonal_expansion_timepoint_plots()
#get_tops_plots(plot_folder = plot_folder)



