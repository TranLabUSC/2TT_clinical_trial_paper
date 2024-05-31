# subset the single cell data for itt patients and remove 293t cells

library(Seurat)
library(ggplot2)

seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/2tt_jci_combined_seurat_object.rds")
# remove IDH positive patients
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")
gbm_patients <- survival_df[survival_df$IDH.1.mutation == "negative", "Subject.ID"]
seurat_object_all_cells <- subset(seurat_object_all_cells, subset = donor %in% gbm_patients)
# remove tumor samples
seurat_object_all_cells <- subset(seurat_object_all_cells, subset = TimePoint != "_")
# remove timepoint 5 patient 19
seurat_object_all_cells <- subset(seurat_object_all_cells, subset = origin != "5_19")
# remove 293t cells
clusters_293t <- c(9,10,13,15)
seurat_clusters_to_keep <- as.numeric(setdiff(unique(seurat_object_all_cells@meta.data$seurat_clusters), clusters_293t))
seurat_object_all_cells <- subset(seurat_object_all_cells, subset = seurat_clusters %in% seurat_clusters_to_keep)

saveRDS(seurat_object_all_cells, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")





seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/2tt_jci_combined_seurat_object_t_cell_res_1_nk_removed.rds")
# remove IDH positive patients
survival_data <- "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/patient_info_2TT_survival.csv"
survival_df <- read.csv(survival_data, sep = "\t")
gbm_patients <- survival_df[survival_df$IDH.1.mutation == "negative", "Subject.ID"]
seurat_object_all_cells <- subset(seurat_object_t_cells, subset = donor %in% gbm_patients)
# remove tumor samples
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = TimePoint != "_")
# remove patient 19 timepoint 5
seurat_object_t_cells <- subset(seurat_object_t_cells, subset = origin != "5_19")

saveRDS(seurat_object_t_cells, "/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_T_Cells_seurat_object_GBM.rds")



seurat_object_all_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/GBM/2TT_seurat_object_GBM.rds")
seurat_object_t_cells <- readRDS("/project/dtran642_927/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate_ge_only/Nextflow/Aggregated/outs/NETZEN_analysis/Finale/ITT/2TT_T_Cells_seurat_object_ITT.rds")
p <- DimPlot(seurat_object_all_cells, label = FALSE, raster = FALSE)
p + coord_flip()