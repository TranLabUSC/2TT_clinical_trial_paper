setwd("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code")
source("get_metadata.R")
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/"
)
#source("code/analysis.R")

library(immunarch)
library(cowplot)




#data_folder = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/ClonotypingData"

#clonotyping_data = repLoad(data_folder)

track_patient_clone = function(patient = "p18",
                               top,
                               plot_folder,
                               metadata,
                               trackClonotype_df,
                               plot_with_per_sample = 250,
                               plot_height_per_plot = 800)
{
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  
  #tc2  =trackClonotypes(clonotyping_data$data[patient_samples$sample_id], list(nrow(patient_samples), 1000), .col = "aa")
  tc1  = trackClonotype_df[, patient_samples$sample_id]
  # Eliminate empty rows:
  tc1$total = rowSums(tc1)
  tc1  = tc1 %>% filter(total > 0)
  tc1 = tc1[, 1:(ncol(tc1) - 1)]
  tc1[, "sequence"] = rownames(tc1)
  tc1 = tc1[, c(ncol(tc1), 1:(ncol(tc1) - 1))]
  
  make_clonal_tracking_plot = function(clonal_tracking_data, plot_label)
  {
    class(clonal_tracking_data) = c("immunr_dynamics", "data.table", "data.frame")
    p = vis(clonal_tracking_data) + theme_bw() +
      theme(
        legend.position = "none" ,
        text = element_text(size = 40),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 40, angle = 90)
      ) +  scale_y_continuous(limits = c(0, 1))  +
      scale_x_discrete(labels =  patient_samples$integrated_sample_type)  +
      ggtitle(plot_label) + 
      theme(plot.title = element_text(hjust = 0.5,size = 40))  
      # annotate(
      #   geom = "text",
      #   x = 1,
      #   y = 0.9,
      #   label = plot_label,
      #   size = 20,
      #   hjust = 0
      # )
    return(p)
  }
  
  make_clonal_tracking_plot_per_sample_id = function(sample_id)
  {
    tc_sorted = tc1[order(tc1[, sample_id], decreasing = TRUE),]
    top_sequences = tc_sorted[1:top, 1]
    integrated_sample_type = as.character(patient_samples[patient_samples$sample_id == sample_id, "integrated_sample_type"])
    class(tc_sorted) = c("immunr_dynamics", "data.table", "data.frame")
    plot_label = paste(patient, "\n", "top", top , "of", integrated_sample_type)
    p = make_clonal_tracking_plot(tc_sorted[1:top, ], plot_label)
    return(list(p, top_sequences))
  }
  
  make_plot_list = function()
  {
    plot_width = plot_with_per_sample * max(6, nrow(patient_samples))
    plot_height = plot_height_per_plot  + max(nchar(as.vector(
      patient_samples$integrated_sample_type
    )))
    plot_list =  list()
    combined_top_sequences = vector()
    for (sample_id in patient_samples$sample_id)
    {
      out = make_clonal_tracking_plot_per_sample_id(sample_id = sample_id)
      p = out[[1]]
      top_sequences = out[[2]]
      integrated_sample_type = as.character(patient_samples[patient_samples$sample_id == sample_id, "integrated_sample_type"])
      print(integrated_sample_type)
      filename = paste0(
        plot_folder,
        "/",
        patient,
        "_top",
        top ,
        "_",
        integrated_sample_type,
        "_",
        sample_id,
        ".png"
      )
      png(
        file =  filename ,
        units = "px",
        width = plot_width,
        height = plot_height
      )
      print(p)
      dev.off()
      
      plot_list[[integrated_sample_type]] = p
      combined_top_sequences = c(combined_top_sequences, top_sequences)
    }
    combined_top_sequences = unique(combined_top_sequences)
    return(list(plot_list = plot_list, combined_top_sequences = combined_top_sequences))
  }
  
  clonal_plots = make_plot_list()
  clonal_plot_list = clonal_plots[["plot_list"]]
  combined_top_sequences = clonal_plots[["combined_top_sequences"]]
  make_combined_top_sequences_plot = function(combined_top_sequences)
  {
    tc_combined = tc1[combined_top_sequences,]
    class(tc_combined) = c("immunr_dynamics", "data.table", "data.frame")
    plot_label =  paste(patient, "\n", "top", top , "of all samples")
    combined_top_sequences_plot = make_clonal_tracking_plot(clonal_tracking_data = tc_combined, plot_label = plot_label)
    filename =  paste0(plot_folder, "/", patient,  "_top", top , "_all_samples.png")
    png(
      file =  filename ,
      units = "px",
      width = plot_with_per_sample * nrow(patient_samples),
      height =  plot_height_per_plot  + max(nchar(
        as.vector(patient_samples$integrated_sample_type)
      ))
    )
    print(combined_top_sequences_plot)
    dev.off()
    
    return(combined_top_sequences_plot)
  }
  combined_top_sequences_plot = make_combined_top_sequences_plot(combined_top_sequences)
  
  layout_all_plots = function(combined_top_sequences_plot,
                              clonal_plot_list)
  {
    filename = paste0(plot_folder, "/", patient,  "_top", top , "_allsamples.png")
    plot_list = append(list(combined_top_sequences_plot), clonal_plot_list)
    combined_plot = plot_grid(plotlist = plot_list,
                              align = "h",
                              ncol = 1)
    # layout multiplot with cowplot
    filename = paste0(plot_folder, "/", patient,  "_top", top , "_combined.png")
    plot_width = plot_with_per_sample * max(6, nrow(patient_samples))
    plot_height =  (nrow(patient_samples) + 1) * (plot_height_per_plot  + max(nchar(
      as.vector(patient_samples$integrated_sample_type)
    )))
    
    png(
      file =  filename ,
      units = "px",
      width = plot_width,
      height = plot_height
    )
    print(combined_plot)
    dev.off()
    return(combined_plot)
  }
  
  
  combined_plot = layout_all_plots(combined_top_sequences_plot, clonal_plot_list)
  
  out = list(
    combined_plot = combined_plot,
    combined_top_sequences_plot = combined_top_sequences_plot,
    clonal_plot_list = clonal_plot_list
  )
  return(out)
}




make_tumor_volume_plot = function(patient = "p18")
{
  print()
  
}


layout_all_plots = function(patient,
                            clonal_tracking_plot,
                            diversity_plot,
                            tumor_volume_plot)
{
  print()
  
}


make_clonal_plot_all_patients = function(plot_folder = "outs",
                                         top = 100,
                                         trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt",
                                         archive_file =  "clonal_plot_all_patients.RDS")
{
  metadata = get_metadata_sorted()
  trackClonotype_df <-
    read.delim(trackClonotype_df_file, check.names = FALSE)
  
  
  clonal_plot_all_patients = list()
  for (patient in unique(metadata$Patient))
  {
    print(paste("creating plots for patient", patient))
    plot_folder_patient = paste0(plot_folder, "/", patient)
    dir.create(plot_folder_patient)
    
    tryCatch({
      clonal_plot = track_patient_clone(
        patient = patient,
        top = top,
        plot_folder = plot_folder_patient,
        metadata,
        trackClonotype_df,
        plot_with_per_sample = 250,
        plot_height_per_plot = 800
      )
      clonal_plot_all_patients[[patient]] = clonal_plot
    }
    
    , error = function(cond)
    {
      print(cond)
    })
  }
  
  saveRDS(clonal_plot_all_patients, file = archive_file)
  
}


make_clonal_plot_all_patients(plot_folder = "outs",
                              top = 10,
                              trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt",
                              archive_file =  "top10_clonal_plot_all_patients.RDS")
make_clonal_plot_all_patients(plot_folder = "outs",
                              top = 20,
                              trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt",
                              archive_file =  "top20_clonal_plot_all_patients.RDS")
make_clonal_plot_all_patients(plot_folder = "outs",
                              top = 100,
                              trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt",
                              archive_file =  "top100_clonal_plot_all_patients.RDS")

