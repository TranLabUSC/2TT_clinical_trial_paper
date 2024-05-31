setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("trackclone_plots.R")
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)
library(gridExtra)
library(cowplot)
library()
prepare_data = function(patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv")
{
  metadata = get_metadata_sorted()
  patient_survival_data <- read.delim(patient_survival_data_file)
  patient_survival_data$Patient = paste0("p", patient_survival_data$Subject.ID)
  patients = unique(metadata$Patient)
  setdiff(patients, patient_survival_data$Patient) # Missing info for patient: 13, 27, 34
  return(list(metadata = metadata, patient_survival_data = patient_survival_data))
}


make_pfs_os_plot = function(patient = "p16",
                            plot_width = 1000,
                            plot_height = 200,
                            metadata,
                            plot_folder,
                            patient_survival_data)
  
{
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  survival_data = patient_survival_data %>% filter(Patient == patient)
  if (nrow(survival_data) == 0)
  {
    return(NULL)
  }
  ORR = survival_data$Best.ORR
  col_dict = list(
    "PD" = "red",
    "SD" = "yellow",
    "PR" = "green",
    "PR/near CR" = "blue"
  )
  fill_color = col_dict[[ORR]]
  
  OS_PFS = survival_data %>% select(PFS.from.Diagnosis..months., OS.from.diagnosis..months.)
  OS_PFS = t(OS_PFS)
  OS_PFS = as.data.frame(OS_PFS)
  OS_PFS$variable = c("PFS", "OS")
  
  
  
  p1 = ggplot(data =  OS_PFS, aes(x = variable, y =
                                    V1, label = V1)) + geom_bar(stat = "identity", fill = fill_color) +
    geom_text(hjust = -0.2 ,
              size = 10) +   labs(y = "months since diagnosis" , x = "") + coord_flip() +
    theme(
      legend.position = "none" ,
      text = element_text(size = 30),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.title = element_blank(),
      axis.text = element_text(size = 40),
      axis.title = element_text(size = 40)
    ) + scale_y_continuous(limits = c(
      0,
      max(patient_survival_data$OS.from.diagnosis..months.) + 5
    ))
  print(p1)
  
  
  filename = paste0(plot_folder, "/", patient, "_OS_PFS.png")
  
  
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


make_clinical_info_table = function(patient = "p16",
                                    plot_with_per_sample = 250,
                                    plot_height_per_plot = 1000,
                                    metadata,
                                    plot_folder,
                                    patient_survival_data)
{
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  survival_data = patient_survival_data %>% filter(Patient == patient)
  if (nrow(survival_data) == 0)
  {
    return(NULL)
  }
  survival_data_t = as.data.frame(t(survival_data))
  select_clinical_data = c(
    "Patient" ,
    "Age",
    "Sex",
    "MGMT.methylation",
    "IDH.1.mutation" ,
    "Best.ORR" ,
    "Tumor.location"  ,
    "Dead"
  )
  survival_data_t = as.data.frame(survival_data_t[select_clinical_data,])
  survival_data_t$clinical.indicator =  c(
    "Patient" ,
    "Age",
    "Sex",
    "MGMT.methylation",
    "IDH.1.mutation" ,
    "Best.ORR" ,
    "Tumor.location"  ,
    "Dead"
  )
  survival_data_t = survival_data_t[, c(2, 1)]
  colnames(survival_data_t) = c("clinical", "value")
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params = list(cex = 4.0)),
    colhead = list(fg_params = list(cex = 4.0)),
    rowhead = list(fg_params = list(cex = 4.0))
  )
  tbl <- tableGrob(survival_data_t, rows = NULL, theme = mytheme)
  
  
  
  
  return(tbl)
  
}


make_plot_for_all_patients = function(plot_folder = "outs",
                                      plot_width = 1000,
                                      plot_height = 800)
{
  prepared_data = prepare_data(patient_survival_data_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv")
  metadata = prepared_data[["metadata"]]
  patient_survival_data = prepared_data[["patient_survival_data"]]
  patients = unique(metadata$Patient)
  patient_info_plots = list()
  for (patient in patients)
  {
    print(paste("creating plots for patient", patient))
    plot_folder_patient = paste0(plot_folder, "/", patient)
    dir.create(plot_folder_patient)
    #tryCatch({
    
    survival_plot = make_pfs_os_plot(
      patient = patient,
      plot_width = 1000,
      plot_height = 200,
      metadata = metadata,
      plot_folder = plot_folder_patient,
      patient_survival_data = patient_survival_data
    )
    clinical_table = make_clinical_info_table(
      patient = patient,
      metadata = metadata,
      plot_folder = plot_folder_patient,
      patient_survival_data = patient_survival_data
    )
    if (!is.null(survival_plot) & !is.null(clinical_table))
    {
      combined_clinical_plot = plot_grid(
        clinical_table,
        survival_plot ,
        nrow = 2,
        ncol = 1,
        rel_widths = c(1, 1),
        rel_heights = c(3, 1)
      )
      filename = paste0(plot_folder_patient,
                        "/",
                        patient,
                        "_combined_clinical.png")
      
      png(
        file =  filename ,
        units = "px",
        width = plot_width,
        height = plot_height
      )
      print(combined_clinical_plot)
      dev.off()
      
      
      patient_info_plots[[patient]] = list(
        survival_plot = survival_plot,
        clinical_table = clinical_table,
        combined_clinical_plot = combined_clinical_plot
      )
      #}
      
      #, error = function(cond)
      #{
      #  print(cond)
      #})
      #}
    }
  }
  saveRDS(patient_info_plots, file = "patient_info_plots_all_patients.RDS")
  
}
make_plot_for_all_patients()
