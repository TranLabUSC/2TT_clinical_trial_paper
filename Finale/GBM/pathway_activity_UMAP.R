library(cowplot)
library(Seurat)
library(Matrix)
library(future)
options(future.globals.maxSize = 64000 * 1024 ^ 2)
library(patchwork)
library(ggplot2)
library(RColorBrewer)


# ====================================================================================================================
seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_object_all_cells <- NormalizeData(seurat_object_all_cells, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE, assay = "RNA")

seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_T_Cells_seurat_object_GBM.rds")


library(cowplot)
plotPathwayUMAP <- function(gmt_file, pathway_name, seurat_obj, survival_data, file_name) {
  # Read GMT file and extract genes for the specified pathway
  gmt <- readLines(gmt_file)
  pathway_genes <- NULL
  for (line in gmt) {
    split_line <- strsplit(line, "\t")[[1]]
    if (split_line[1] == pathway_name) {
      pathway_genes <- split_line[-c(1,2)]  # Assuming the first two columns are pathway name and description
      break
    }
  }
  
  # Check if pathway was found
  if (is.null(pathway_genes)) {
    stop("Pathway not found in GMT file.")
  }
  
  # Read survival data and sort patients by survival
  survival_df <- read.csv(survival_data, sep = "\t")
  survival_df <- survival_df[survival_df$IDH.1.mutation == "negative",]
  sorted_patients <- survival_df[order(-survival_df$OS.from.enrollment..months.), 'Subject.ID']
  print(sorted_patients)
  
  # Create a grid of feature plots
  plot_list <- list()
  timepoints <- c("N", "Y", "1", "2", "3", "4", "5", "R", "R2")
  plot_counter <- 1
  for (patient in sorted_patients) {
    print(patient)
    patient_specific_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data$donor == patient,])
    if (length(patient_specific_cells) == 0) {
      next
    }
    patient_seurat_obj <- subset(seurat_obj, cells = patient_specific_cells)
    # Get the genes present in the Seurat object
    seurat_genes <- rownames(patient_seurat_obj@assays$RNA@counts)
    
    # Find intersection of pathway genes and Seurat object genes
    common_genes <- intersect(pathway_genes, seurat_genes)
    
    # Check if there are any common genes
    if (length(common_genes) == 0) {
      stop("None of the pathway genes are found in the Seurat object.")
    }
    
    # Calculate average expression of pathway genes per cell
    # Note: This assumes that the data is already normalized
    pathway_avg_expression <- apply(GetAssayData(patient_seurat_obj, assay = "RNA", slot = "data")[common_genes, ], 2, mean, na.rm = TRUE)
    
    # Check for cells at timepoint N
    timepoint_N_barcodes <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == "N", ])
    
    if (length(timepoint_N_barcodes) > 0) {
      # Calculate mean and standard deviation for cells at timepoint N
      mean_N <- mean(pathway_avg_expression[timepoint_N_barcodes])
      sd_N <- sd(pathway_avg_expression[timepoint_N_barcodes])
      
      # Scale based on timepoint N
      scaled_pathway_avg_expression <- (pathway_avg_expression - mean_N) / sd_N
    } else {
      # Use scale function if no cells at timepoint N
      scaled_pathway_avg_expression <- scale(pathway_avg_expression)[, 1]
    }
    
    
    pathway_avg_expression <- scaled_pathway_avg_expression
    names(pathway_avg_expression) <- colnames(patient_seurat_obj)
    
    # Add this as a metadata column
    patient_seurat_obj[["pathway_avg_expression"]] <- pathway_avg_expression
    
    # Find the minimum number of cells across all timepoints for this patient, excluding timepoints with no cells
    cell_counts <- sapply(timepoints, function(tp) {
      sum(patient_seurat_obj$TimePoint == tp)
    })
    print(cell_counts)
    min_cells <- min(cell_counts[cell_counts > 0])  # Exclude timepoints with no cells
    
    for (tp in timepoints) {
      patient_tp_cells <- rownames(patient_seurat_obj@meta.data[patient_seurat_obj@meta.data$TimePoint == tp, ])
      
      if (length(patient_tp_cells) > 0) {
        # Randomly sample cells if this timepoint has more cells than min_cells
        # if (length(patient_tp_cells) > min_cells) {
        #   set.seed(123)  # For reproducibility
        #   patient_tp_cells <- sample(patient_tp_cells, min_cells)
        # }
        print(length(patient_tp_cells))
        # patient_seurat_obj <- subset(seurat_obj, cells = patient_cells)
        p <- FeaturePlot(patient_seurat_obj, features = "pathway_avg_expression", cells = patient_tp_cells) &
          scale_colour_gradientn(limits = c(-4, 4), colours = rev(brewer.pal(n = 11, name = "RdBu"))) &
          theme(legend.position = 'none',
                plot.title = element_blank(), 
                axis.title = element_blank(), 
                panel.border = element_rect(colour = "black", fill=NA, size=1), 
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) # Adjust the numbers as needed) + # This adds the boxed appearance and removes titles
        plot_list[[plot_counter]] <- p + coord_flip()
      } else {
        plot_list[[plot_counter]] <- ggplot() + 
          geom_blank() + 
          scale_x_continuous(limits = c(-5, 5), breaks = c(-5, 0, 5)) +
          scale_y_continuous(limits = c(-10, 10), breaks = c(-10, -5, 0, 5, 10)) +
          theme_minimal() +
          theme(legend.position = 'none',
                plot.title = element_blank(), 
                axis.title = element_blank(), 
                axis.text = element_text(size = 12, colour = "black"), # Adjust the axis text size
                axis.ticks = element_line(colour = "black"), # Ensure ticks are visible
                panel.border = element_rect(colour = "black", fill=NA, size=1), 
                panel.grid.major = element_blank(), # Remove major grid lines
                panel.grid.minor = element_blank(), # Remove minor grid lines
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")) 
      }
      plot_counter <- plot_counter + 1
    }
  }
  
  # Combine plots into a grid using cowplot
  plot_grid <- plot_grid(plotlist = plot_list, ncol = length(timepoints), align = 'v', rel_widths = c(0.5), rel_heights = c(0.5))
  
  # Save the grid plot as a PDF
  ggsave(filename = file_name, plot = plot_grid, device = "pdf", width = 4 * length(timepoints), height = 4 * length(sorted_patients), limitsize = FALSE)
}

# Example usage

survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
gmt_file <- "/project/dtran642_927/SonLe/USC_Source/source/NetZenPathwayAnalysis/subsubset_output.gmt"
pathways <- c("GO_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS", "GO_RESPONSE_TO_TYPE_I_INTERFERON")
# pathways <- c("GO_ADAPTIVE_IMMUNE_RESPONSE", "GO_T_CELL_ACTIVATION", "GO_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE", "GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS", "REACTOME_IMMUNE_SYSTEM", "GO_IMMUNE_RESPONSE-ACTIVATING_SIGNALING_PATHWAY", "GO_ACTIVATION_OF_IMMUNE_RESPONSE", "GO_REGULATION_OF_T_CELL_ACTIVATION", "GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")
for (pathway_name in pathways){
  print(pathway_name)
  file_name <- paste0("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/output_file_", pathway_name, "_all_cells.pdf")
  plotPathwayUMAP(gmt_file, pathway_name, seurat_object_all_cells, survival_data, file_name)
}