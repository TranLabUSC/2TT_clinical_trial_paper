#install.packages("dplyr")
library(dplyr)
# Annotate clones from sc manualy from clonotypes.csv
# Make dictionary mapping TRB to clonotype
# Read bulk clonotyping data by immunoarch
# Relate bulk clonotyping to sc clone by TRB look up, if the TRB does not exist, make new clone.


# Annotate clones from sc manualy from clonotypes.csv

#filtered_contig_annotations <- read.csv("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/ClonotypingData/sc_aggregate_out/filtered_contig_annotations.csv")
#clonotypes <- read.csv("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/ClonotypingData/sc_aggregate_out/clonotypes.csv")




get_clone_dict = function(clonotypes, col = "aa")
{
  real_col = paste0("cdr3s_", col)
  clone_dict = list()
  for (i in 1: nrow(clonotypes))
  {
    clonotype = clonotypes[i, "clonotype_id"]
    cdr3 = as.character(clonotypes[i,real_col])
    chains = unlist(strsplit(cdr3,split=";"))
    for (chain in chains)
    {
      chain_parsed = unlist(strsplit(chain, ":"))
      if (chain_parsed[1]=="TRB")
      {
        sequence = chain_parsed[2]
        if (is.null(clone_dict[[sequence]])) # double bracket
        {
          clone_dict[[sequence]] = as.vector(clonotype) # double bracket to update element
          
        }else
          
        {
          updated_value = c( clone_dict[[sequence]], clonotype)
          clone_dict[[sequence]] = updated_value
          print("clone already exist")
          print( clone_dict[[sequence]])
          if (length(clone_dict[[sequence]]) > 5)
          {
            stop()
          }

        }
      }
      
        
    }
    
  }
  return(clone_dict)
  
}
  
clone_dict = get_clone_dict(clonotypes = clonotypes)



# Make dictionary mapping TRB to clonotype
# Read bulk clonotyping data by immunoarch
# Relate bulk clonotyping to sc clone by TRB look up, if the TRB does not exist, make new clone.




setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/2TT_aggregate/outs/vdj_t"
)

get_cell_fraction_from_scRNAseq = function(filtered_contig_annotations_file = "filtered_contig_annotations.csv")
{

  
  filtered_contig_annotations <-
    read.csv(filtered_contig_annotations_file, stringsAsFactors = FALSE)
  patients = unique(filtered_contig_annotations$Meta.Patient)
  
  patient_clonotypes_with_distict_barcode = filtered_contig_annotations %>% select(c (
    Meta.Patient,
    Meta.Sample.Type,
    origin,
    barcode,
    raw_clonotype_id
  )) %>% distinct()
  patient_clonotypes_number_cell_per_clone_for_each_sample = patient_clonotypes_with_distict_barcode %>% group_by(Meta.Patient, Meta.Sample.Type, origin, raw_clonotype_id) %>% summarise(number_cells = n())
  patient_clonotypes_per_sample_fraction = patient_clonotypes_number_cell_per_clone_for_each_sample %>% group_by(Meta.Patient, Meta.Sample.Type, origin) %>% mutate(
    total_number_cells_in_sample = sum(number_cells),
    clonal_fraction = number_cells / sum(number_cells)
  )
  sample_clonal_fraction_sorted = patient_clonotypes_per_sample_fraction[with(
    patient_clonotypes_per_sample_fraction,
    order(Meta.Patient, Meta.Sample.Type, -clonal_fraction)
  ), ]
  
  sample_clonal_fraction_sorted$time_point = lapply(sample_clonal_fraction_sorted$origin, function(x) {
    strsplit(x, "_")[[1]][1]
  })
  unique(sample_clonal_fraction_sorted$time_point)
  sample_clonal_fraction_sorted$time_point = factor(
    sample_clonal_fraction_sorted$time_point,
    levels = c ("T", "N", "Y", "1", "2", "3", "4", "5", "R", "R2")
  )
  sample_clonal_fraction_sorted = sample_clonal_fraction_sorted[with(
    sample_clonal_fraction_sorted,
    order(Meta.Patient, time_point, -clonal_fraction)
  ), ]
  write.table(
    sample_clonal_fraction_sorted,
    file = "sample_clonal_fraction_sorted.txt",
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  #sample_clonal_fraction_sorted = sample_clonal_fraction_sorted[, c("Meta.Patient", "Meta.Sample.Type", "time_point", "origin", "clonal_fraction", "number_cells", "total_number_cells_in_sample")]
  
  
}