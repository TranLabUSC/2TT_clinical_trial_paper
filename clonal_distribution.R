setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("get_metadata.R")
library(ggplot2)
library(cowplot)
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)




make_plot = function(Clonotype_df = Clonotype_df, sample_name = "38T0_sc", top = 50 , plot_label = "clonal distribution" )
{
print(sample_name)
clones = Clonotype_df[, sample_name]

clones = sort(clones, decreasing = TRUE)
top_clones = clones[1:top]
top_clones = data.frame(clone = 1:top, proportion = top_clones)
p1 = ggplot(data =  top_clones) + geom_col(aes(x = clone, y =
                                                      proportion)) +
  theme_bw() +
  theme(
    legend.position = "none" ,
    text = element_text(size = 30),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +  
  scale_y_continuous(limits =c(0, 0.3))  +
  ggtitle(plot_label) +
  theme(plot.title = element_text(
    hjust = 0.5,
    vjust = 0.2,
    size = 40
  ))
print(p1)
return(p1)
}

make_patient_clonal_distribution_plot = function(Clonotype_df = Clonotype_df,
                                                 patient = "p9",
                                                 top = 50 ,
                                                 plot_folder = "outs",
                                                 metadata = metadata)
{
  
  
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  plot_list = list()
  patient_samples$sample_id = factor(patient_samples$sample_id, levels =  patient_samples$sample_id)
  for (sample_ in sample_order)
  {
    row = patient_samples[patient_samples$integrated_sample_type == sample_,]
    if (nrow(row) == 0)
    {
      plot_list[sample_] = list(NULL)
    }else
    {
    sample_id = as.character(row[,"sample_id"][1])
    plot_label = paste( patient, "clonal_distribution", "\n", sample_)
    p = make_plot(Clonotype_df = Clonotype_df, 
                  sample_name = sample_id, 
                  top = top,
                  plot_label = plot_label )
    plot_list[[sample_]] = p
    }
  
    
    
  }
  
  combined_plot = plot_grid(plotlist=plot_list, ncol =1 )
 
 filename = paste0(plot_folder, "/", patient, "_clonal_distribution.pdf")
 plot_width = 15
 plot_height =  5 * length(plot_list)
 cairo_pdf(
   file =  filename ,
   width = plot_width,
   height = plot_height
 )
 print(combined_plot)
 dev.off()
 #return(p1)
 saveRDS(list(plot_list = plot_list, combined_plot = combined_plot), paste0(plot_folder, "/", patient, "_clonal_distribution.RDS"))
return(combined_plot)
   
}



make_plot_for_all_patients = function(plot_folder = "outs", metadata = metadata, top = 50)
{
  Clonotype_df <- read.delim("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt", 
                           check.names = FALSE)
  metadata = get_metadata_sorted()
  patients = unique(metadata$Patient)
  clonal_dist_plots = list()
  for (patient in patients)
  {
    print(paste("creating plots for patient", patient))
    plot_folder_patient = paste0(plot_folder, "/", patient)
    dir.create(plot_folder_patient)
    combined_plot =  make_patient_clonal_distribution_plot(Clonotype_df = Clonotype_df,
                                                                      patient = patient,
                                                                      top = top,
                                                                      plot_folder = plot_folder_patient,
                                                           metadata = metadata)
   clonal_dist_plots[[patient]] = combined_plot 
  }
  saveRDS(clonal_dist_plots, file = "clonal_distribution_plots_all_patients.RDS")
}

make_plot_for_all_patients()

