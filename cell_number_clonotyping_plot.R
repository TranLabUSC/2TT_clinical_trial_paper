setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("trackclone_plots.R")
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)


library(cowplot)

prepare_data = function()
  
{
metadata = get_metadata_sorted()
filtered_contig_annotations <-
  read.csv(
    "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/sc_aggregate_out/filtered_contig_annotations.csv"
  )
trackClonotype_df.abs <-
  read.delim(
    "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.abs.txt",
    check.names = FALSE
  )
number_cells_per_sample = as.data.frame(colSums(trackClonotype_df.abs[, 2:ncol(trackClonotype_df.abs)]))
colnames(number_cells_per_sample) = "number_cell"
number_cells_per_sample$sample = rownames(number_cells_per_sample)

number_clones_per_sample  = trackClonotype_df.abs[, 2:ncol(trackClonotype_df.abs)]
number_clones_per_sample[number_clones_per_sample  > 0] = 1
number_clones_per_sample = colSums(number_clones_per_sample)
number_clones_per_sample = data.frame(sample = names(number_clones_per_sample),
                                      number_clones = number_clones_per_sample)
number_cells_clones_per_sample = merge(
  number_cells_per_sample,
  number_clones_per_sample,
  by = "sample",
  all = TRUE,
  sort = FALSE
)
row.names(number_cells_clones_per_sample) = number_cells_clones_per_sample$sample


number_clones_with_multi_cells  = trackClonotype_df.abs[, 2:ncol(trackClonotype_df.abs)]
number_clones_with_multi_cells[number_clones_with_multi_cells  < 2] = 0
number_clones_with_multi_cells[number_clones_with_multi_cells  >= 2] = 1
number_clones_with_multi_cells = colSums(number_clones_with_multi_cells)
number_clones_with_multi_cells = data.frame(
  sample = names(number_clones_with_multi_cells),
  number_multi_cell_clones = number_clones_with_multi_cells
)

number_cells_clones_per_sample = merge(
  number_cells_per_sample,
  number_clones_per_sample,
  by = "sample",
  all = TRUE,
  sort = FALSE
)
number_cells_clones_per_sample = merge(
  number_cells_clones_per_sample,
  number_clones_with_multi_cells,
  by = "sample",
  all = TRUE,
  sort = FALSE
)
row.names(number_cells_clones_per_sample) = number_cells_clones_per_sample$sample
number_cells_clones_per_sample$number_single_cell_clones = number_cells_clones_per_sample$number_clones - number_cells_clones_per_sample$number_multi_cell_clones


number_cells_clones_per_sample$cells_per_clone = number_cells_clones_per_sample$number_cell /
  number_cells_clones_per_sample$number_clones
write.table(
  number_cells_clones_per_sample,
  "number_cells_clones_per_sample.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

return(list(metadata= metadata, all_patients_data = number_cells_clones_per_sample))

}

# Make the plot of cell numbers or clone number:
make_plot = function(patient = "p19",
                     plot_width = 1000,
                     plot_height = 200,
                     metadata,
                     plot_folder,
                     all_patients_data,
                     color_cutoff = 500,
                     y_limit = 10000,
                     x = "sample",
                     y = "number_clones")
  # data type can be "clone" or "cell"
{
  patient_samples = metadata %>% filter(Patient == patient, dataset != "bulk")
  patient_data = all_patients_data[patient_samples$sample_id,]
  patient_data$sample = factor(patient_data$sample, levels = patient_data$sample)
  patient_data$fill_color = sapply(patient_data[, y], function(x) {
    if (x >= color_cutoff)
      return("blue")
    else
    {
      return("red")
    }
  })
  
  title_dict = list(
    number_cell = " number of cells",
    number_clones  = "number of clones",
    number_multi_cell_clones = "number of multi-cell clones",
    number_single_cell_clones = "number of single-cell clones"
  )
  plot_label = paste0(patient, " - ", title_dict[[y]])
  p1 = ggplot(data = patient_data) + geom_col(aes(x = .data[[x]], y =
                                                    .data[[y]]), fill = patient_data$fill_color) +
    geom_text(
      aes(x = .data[[x]], y = .data[[y]], label = .data[[y]]),
      hjust  = -0.1,
      size = 10,
      angle = 90
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
      axis.text.x = element_text(
        size = 40,
        angle = 90,
        color = patient_data$fill_color
      )
    ) +  scale_y_continuous(limits = c(0, y_limit))  +
    scale_x_discrete(labels =  patient_samples$integrated_sample_type)  +
    ggtitle(plot_label) +
    theme(plot.title = element_text(hjust = 0.5, size = 40))
  print(p1)
  
  
  filename = paste0(plot_folder, "/", patient, "_", y, ".png")
  
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

make_patient_combined_plot = function(patient,
                             metadata,
                             all_patients_data,
                             plot_folder = ".",
                             plot_width = 1000,
                             plot_height = 800)
{
  cell_number_plot = make_plot(
    patient = patient,
    plot_width = 1000,
    plot_height = 800,
    metadata = metadata,
    plot_folder = plot_folder,
    all_patients_data = all_patients_data ,
    x = "sample",
    y = "number_cell",
    color_cutoff = 500,
    y_limit = 10000
  )
  
  
  clone_number_plot  = make_plot(
    patient = patient,
    plot_width = plot_width,
    plot_height = plot_height,
    metadata = metadata,
    plot_folder = plot_folder,
    all_patients_data = all_patients_data ,
    x = "sample",
    y = "number_clones",
    color_cutoff = 50,
    y_limit = 4000
  )
  number_multi_cell_clones_plot  = make_plot(
    patient = patient,
    plot_width = plot_width,
    plot_height = plot_height,
    metadata = metadata,
    plot_folder = plot_folder,
    all_patients_data = all_patients_data ,
    x = "sample",
    y = "number_multi_cell_clones",
    color_cutoff = 10,
    y_limit = 4000
  )
  
  number_single_cell_clones_plot  = make_plot(
    patient = patient,
    plot_width = plot_width,
    plot_height = plot_height,
    metadata = metadata,
    plot_folder = plot_folder,
    all_patients_data = all_patients_data ,
    x = "sample",
    y = "number_single_cell_clones",
    color_cutoff = 50,
    y_limit = 4000
  )
  
  combined_plot = plot_grid(
    plotlist = list(
      cell_number_plot,
      clone_number_plot,
      number_multi_cell_clones_plot,
      number_single_cell_clones_plot
    ),
    nrow = 4,
    ncol = 1
  )
  
  filename = paste0(plot_folder, "/", patient, "_sample_cell_numbers.png")
  png(file =  filename ,
      width = plot_width,
      height = plot_height * 4)
  print(combined_plot)
  dev.off()
  
  return(list(combined_plot= combined_plot, 
              cell_number_plot = cell_number_plot, 
              clone_number_plot = clone_number_plot,  
              number_multi_cell_clones_plot = number_multi_cell_clones_plot,
              number_single_cell_clones_plot = number_single_cell_clones_plot))
}

make_plot_for_all_patients = function(plot_folder = "outs",
                                      plot_width = 1000,
                                      plot_height = 800)
{

all_patient_plot_list = list()
prepared_data = prepare_data()

metadata = prepared_data[["metadata"]]
all_patients_data = prepared_data[["all_patients_data"]]
patients = unique(metadata$Patient)
patient_info_plots = list()
for (patient in patients)
{

  print(paste("creating plots for patient", patient))
  plot_folder_patient = paste0(plot_folder, "/", patient)
  dir.create(plot_folder_patient)
 plotlist = make_patient_combined_plot(patient = patient, metadata = metadata, all_patients_data = all_patients_data, plot_folder = plot_folder_patient)
  all_patient_plot_list[[patient]] = plotlist
}

saveRDS(all_patient_plot_list, file = "cell_clone_numbers_plots.RDS")
return(all_patient_plot_list)
}

make_plot_for_all_patients()





