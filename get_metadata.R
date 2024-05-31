library(dplyr)
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


get_metadata_sorted = function(meta_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/Metadata/combined_meta.txt")
{
  combined_meta <- read.delim(meta_file)
  combined_meta$integrated_sample_type = paste(combined_meta$Sample.Type,
                                               combined_meta$timepoint,
                                               combined_meta$dataset,
                                               sep = "_")
  combined_meta = combined_meta %>% filter(dataset != "bulk")
  combined_meta$integrated_sample_type = factor(combined_meta$integrated_sample_type, levels = sample_order)
  metadata = combined_meta %>% arrange(Patient, integrated_sample_type)
  
  return(metadata)
}

to_bool = function(x)
{
  if (x == "Y")
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}