setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("get_metadata.R")
library(ggplot2)
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)


prepare_data = function(tumor_size_file =  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/tumor_size_2TT.csv" , sample_order = sample_order)
{
  metadata = get_metadata_sorted()
  metadata$original_sample_id = sapply(metadata$sample_id, function(x)
    unlist(strsplit(x, split = "_"))[[1]])
  
  patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
  patient_survival_data <- read.delim(patient_survival_data_file)
  patient_survival_data$Patient = paste0("p", patient_survival_data$Subject.ID)
  patients = unique(metadata$Patient)
  setdiff(patients, patient_survival_data$Patient) # Missing info for patient: 13, 27, 34
  
  tumor_size_df <-
    read.delim(tumor_size_file)
  tumor_size_df$Patient = paste0("p", tumor_size_df$Sample.ID)
  # Tumor size at the time point of corresponding samples
  tumor_size_blood_sample = tumor_size_df[, c("Blood.Sample.Name", "Tumor.Size..cm2.")]
  tumor_size_tumor_sample = tumor_size_df[, c("Tumor.sample.name" , "Tumor.Size..cm2.")]
  colnames(tumor_size_blood_sample) = c("Sample.ID", "Tumor.Size")
  colnames(tumor_size_tumor_sample) = c("Sample.ID", "Tumor.Size")
  tumor_size_combined = rbind(tumor_size_blood_sample, tumor_size_tumor_sample)
  tumor_size_combined = tumor_size_combined %>% filter(Sample.ID != "")
  
  metadata = merge(
    metadata,
    tumor_size_combined,
    by.x = "original_sample_id",
    by.y = "Sample.ID",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )

  
  metadata$integrated_sample_type = factor(metadata$integrated_sample_type, levels = sample_order)
  metadata = metadata %>% arrange(Patient, integrated_sample_type)
  
  
  write.table(
    metadata,
    "metadata_with_tumor_volume.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  return(metadata)
}

make_tumor_size_plot = function(patient = "p16",
                                plot_with_per_sample = 250,
                                plot_height_per_plot = 1000,
                                metadata,
                                plot_folder)
{
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  
  patient_samples$sample_id = factor(patient_samples$sample_id, levels =  patient_samples$sample_id)
  patient_samples$Tumor.Size.NoNA = patient_samples$Tumor.Size
  patient_samples$Tumor.Size.NoNA[is.na(patient_samples$Tumor.Size.NoNA)] = 0
  plot_label = paste(patient, "Tumor size (cm2)")
  p1 = ggplot(data =  patient_samples) + geom_col(aes(x = sample_id, y =
                                                        Tumor.Size)) +
    geom_text(
      aes(x = sample_id, y = Tumor.Size.NoNA, label = Tumor.Size),
      vjust = -0.2,
      size = 10
    ) +
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
    ) +  scale_y_continuous()  +
    scale_x_discrete(labels =  patient_samples$integrated_sample_type)  +
    ggtitle(plot_label) +
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.2,
      size = 40
    ))
  
  filename = paste0(plot_folder, "/", patient, "_tumor_size.png")
  plot_width = plot_with_per_sample * max(6, nrow(patient_samples))
  plot_height =  plot_height_per_plot  + max(nchar(as.vector(
    patient_samples$integrated_sample_type
  )))
  
  png(
    file =  filename ,
    units = "px",
    width = plot_width,
    height = plot_height
  )
  print(p1)
  dev.off()
  return(p1)
  
}



make_plot_for_all_patients = function(plot_folder = "outs",
                                      plot_with_per_sample = 250,
                                      plot_height_per_plot = 1000)
{
  metadata = prepare_data(tumor_size_file =  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/tumor_size_2TT.csv")
  patients = unique(metadata$Patient)
  tumor_size_plots = list()
  for (patient in patients)
  {
    print(paste("creating plots for patient", patient))
    plot_folder_patient = paste0(plot_folder, "/", patient)
    dir.create(plot_folder_patient)
    tryCatch({
      patient_tumor_size_plot = make_tumor_size_plot(
        patient = patient,
        plot_with_per_sample = plot_with_per_sample ,
        plot_height_per_plot = plot_height_per_plot,
        plot_folder = plot_folder_patient,
        metadata = metadata
        
      )
      tumor_size_plots[[patient]] = patient_tumor_size_plot
    }
    
    , error = function(cond)
    {
      print(cond)
    })
  }
  saveRDS(tumor_size_plots, file = "tumor_size_plots_all_patients.RDS")
}

make_plot_for_all_patients()


