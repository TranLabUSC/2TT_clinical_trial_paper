setwd("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code")

source("get_metadata.R")
library(ggplot2)

setwd("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping")

#Prepare data
metadata = get_metadata_sorted()
prepare_data = function(clonal_diversity_file= "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/clonal_diversity.txt")
{

clonal_diversity <- read.delim(clonal_diversity_file)
clonal_diversity$sample_id = rownames(clonal_diversity)
clonal_diversity = clonal_diversity[,c(2,1)] 
clonal_diversity_merged  = merge(metadata, clonal_diversity, by = "sample_id")
clonal_diversity_merged$inversion_diversity 
# inversion_diversity: inversion of Simpson Diversity Index, as a measure of clonal expansion.
clonal_diversity_merged$inversion_diversity  = 1/ clonal_diversity_merged$clonotype_diversity
clonal_diversity_merged$inversion_diversity = round(clonal_diversity_merged$inversion_diversity, 2 )
clonal_diversity_merged$clonotype_diversity = round(clonal_diversity_merged$clonotype_diversity, 2 )
clonal_diversity_merged$sample_id = factor(clonal_diversity_merged$sample_id, levels = clonal_diversity_merged$sample_id )
rownames(clonal_diversity_merged)= clonal_diversity_merged$sample_id
return(clonal_diversity_merged)

}

make_diversity_plot = function(patient = "p16", clonal_diversity_merged, plot_folder,
                               plot_with_per_sample = 250,
                               plot_height_per_plot = 1000, diversity_index = "Simpson")
{
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  clonal_diversity_patient  =  clonal_diversity_merged[patient_samples$sample_id, ]
  clonal_diversity_patient$sample_id = factor(clonal_diversity_patient$sample_id, levels = clonal_diversity_patient$sample_id)
  plot_label = paste(patient, "clonal expansion (inverse of", diversity_index, "diversity index)")
  limit_dict = list(Simpson = c(1, 1.5), Shannon = c(10,0.5))
   p1 = ggplot(data = clonal_diversity_patient) + geom_col( aes(x=sample_id, y=inversion_diversity))+ 
     geom_text(aes(x=sample_id, y=inversion_diversity, label=inversion_diversity), vjust = -0.2, size = 10) +
     geom_line( aes(x= as.numeric(sample_id), y=inversion_diversity) ,stat="identity",color="red",size=2)+ 
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
      axis.text.x = element_text(size = 40, angle = 90)
      ) +  scale_y_continuous(limits = c(0, limit_dict[[diversity_index]][2]))  +
    scale_x_discrete(labels =  patient_samples$integrated_sample_type)  +
     ggtitle(plot_label) + 
     theme(plot.title = element_text(hjust = 0.5,size = 40))  

   plot_label = paste(patient, diversity_index, " Diversity Index")
   p2 = ggplot(data = clonal_diversity_patient) + geom_col( aes(x=sample_id, y=clonotype_diversity))+ 
     geom_text(aes(x=sample_id, y=clonotype_diversity, label=clonotype_diversity), vjust = -0.2, size = 10) +
     geom_line( aes(x= as.numeric(sample_id), y=clonotype_diversity) ,stat="identity",color="blue",size=2)+ 
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
       axis.text.x = element_text(size = 40, angle = 90)
     ) +  scale_y_continuous(limits = c(0, limit_dict[[diversity_index]][1]))  +
     scale_x_discrete(labels =  patient_samples$integrated_sample_type)  +
     ggtitle(plot_label) + 
     theme(plot.title = element_text(hjust = 0.5,size = 40))  
   
   filename = paste0(plot_folder, "/", patient, "_", diversity_index, "_clonal_expansion.png")
   plot_width = plot_with_per_sample * max(6, nrow(patient_samples))
   plot_height =  plot_height_per_plot  + max(nchar(as.vector(patient_samples$integrated_sample_type)))
   
   png(
     file =  filename ,
     units = "px",
     width = plot_width,
     height = plot_height
   )
   print(p1)
   dev.off()
   
   filename = paste0(plot_folder, "/", patient, "_", diversity_index, "_clonal_diversity.png")
   png(
     file =  filename ,
     units = "px",
     width = plot_width,
     height = plot_height
   )
   print(p2)
   dev.off()
   

  return(list(clonal_expansion_plot = p1, clonal_diversity_plot = p2))
}

make_plot_for_all_patients = function(plot_folder = "outs",
                                      clonal_diversity_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/clonal_diversity.txt",
                                      plot_with_per_sample = 250,
                                      plot_height_per_plot = 1000,
                                      diversity_index = "Simpson")
{
  clonal_diversity_merged = prepare_data(clonal_diversity_file=clonal_diversity_file)
  patients = unique(clonal_diversity_merged$Patient)
  diversity_plots = list()
  for (patient in patients)
  {
    print(paste("creating plots for patient", patient))
    plot_folder_patient = paste0(plot_folder, "/", patient)
    dir.create(plot_folder_patient)
    tryCatch({
      patient_diversity_plots = make_diversity_plot(
        patient = patient,
        clonal_diversity_merged = clonal_diversity_merged, 
        plot_folder = plot_folder_patient,
        plot_with_per_sample = plot_with_per_sample ,
        plot_height_per_plot = plot_height_per_plot,
        diversity_index = diversity_index
      )
      diversity_plots[[patient]] = patient_diversity_plots
    }
    
    , error = function(cond)
    {
      print(cond)
    })
  }
  saveRDS(diversity_plots, file = paste0(diversity_index, "_diversity_plots_all_patients.RDS"))
}


make_plot_for_all_patients(clonal_diversity_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/shannon_clonal_diversity.txt",
                           diversity_index = "Shannon")

make_plot_for_all_patients(clonal_diversity_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/simpson_clonal_diversity.txt",
                           diversity_index = "Simpson")
  
