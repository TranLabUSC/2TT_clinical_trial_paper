setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("get_metadata.R")


get_top_sample_clones_change_clonal_timepoints = function( compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                                           metadata,
                                                           tops  = c(10, 20, 50, 100),
                                                            trackClonotype_df, outdir= "clonal_replacement",
                                                           method = "division")
  #Method: subtraction:  clonal_replacement = top_clone_sample2_proportion_at_t2 - top_clone_sample1_proportion_at_t2 
  #       division: clonal_replacement = top_clone_sample2_proportion_at_t2 / max(top_clone_sample1_proportion_at_t2, 0.001) 
{
  outdir = paste0(outdir, "_", method) 
  dir.create(outdir)
  patients = unique(metadata$Patient)
  get_top_sample_clones_change_clonal_replacement_per_patient = function(patient, tops)
  {
    sample1 = metadata %>% filter(Patient == patient,
                                  integrated_sample_type == compared_pair[[1]]) %>% select(sample_id)
    sample1 = sample1$sample_id
    sample2 = metadata %>% filter(Patient == patient,
                                  integrated_sample_type == compared_pair[[2]]) %>% select(sample_id)
    sample2 = sample2$sample_id
    if (length(sample1) == 0 | length(sample2) == 0)
    {
      return(NULL)
    }
    tc1  = trackClonotype_df[, c(sample1, sample2)]
    
    
    get_top_sample_clones_change = function(sample_name, top)
    {
      tc1_sorted_by_sample = tc1[order(tc1[, sample_name], decreasing = TRUE), ]
      top_sample_clones_proportion = colSums(tc1_sorted_by_sample[1:top,])
      top_sample_clones_change = top_sample_clones_proportion[[2]] / top_sample_clones_proportion[[1]]
      return(
        list(
          top_sample_clones_change = top_sample_clones_change ,
          top_sample_clones_proportion = top_sample_clones_proportion
        )
      )
    }
    
    top_sample_clones_change_sample1_list = list()
    top_sample_clones_change_sample2_list = list()
    clonal_replacement_list = list()
    
    for (top in tops)
    {
      top_sample_clones_change_sample1 = get_top_sample_clones_change(sample1, top)
      top_sample_clones_change_sample2 = get_top_sample_clones_change(sample2, top)
      #clonal switch as ratio of top clone changes
      #clonal_replacement = top_sample_clones_change_sample2 / top_sample_clones_change_sample1
      # clonal switch as ratio of top_clone_sample1 / top_clone_sample2 at t2
      top_clone_sample1_proportion_at_t2 = top_sample_clones_change_sample1[["top_sample_clones_proportion"]][[2]]
      top_clone_sample2_proportion_at_t2 = top_sample_clones_change_sample2[["top_sample_clones_proportion"]][[2]]
      if (method == "division")
      {  
      clonal_replacement = top_clone_sample2_proportion_at_t2 / max(top_clone_sample1_proportion_at_t2, 0.001) # add small ammount to avoid inf
      }else
      {
        clonal_replacement = top_clone_sample2_proportion_at_t2 - top_clone_sample1_proportion_at_t2 
      }
      top_sample_clones_change_sample1_list[[paste0("top", top)]] = top_sample_clones_change_sample1[["top_sample_clones_change"]]
      top_sample_clones_change_sample2_list[[paste0("top", top)]] = top_sample_clones_change_sample2[["top_sample_clones_change"]]
      clonal_replacement_list[[paste0("top", top)]] = clonal_replacement
      
    }
    
    return(
      list(
        top_sample_clones_change_sample1_list = top_sample_clones_change_sample1_list,
        top_sample_clones_change_sample2_list = top_sample_clones_change_sample2_list,
        clonal_replacement_list = clonal_replacement_list
      )
    )
    
  }
  # Prepare output dataframes
  top_sample_clones_change_sample1_df = data.frame(matrix(ncol = length(tops), nrow = length(patients)))
  top_sample_clones_change_sample2_df = data.frame(matrix(ncol = length(tops), nrow = length(patients)))
  clonal_replacement_df = data.frame(matrix(ncol = length(tops), nrow = length(patients)))
  rownames(top_sample_clones_change_sample1_df) = patients
  colnames(top_sample_clones_change_sample1_df) = paste0("top_", tops)
  rownames(top_sample_clones_change_sample2_df) = patients
  colnames(top_sample_clones_change_sample2_df) = paste0("top_", tops)
  rownames(clonal_replacement_df) = patients
  colnames(clonal_replacement_df) = paste0("top_", tops)
  
  for (i in 1:length(patients))
  {
    patient = patients[i]
    out = get_top_sample_clones_change_clonal_replacement_per_patient(patient = patient, tops = tops)
    if (!is.null(out))
    {
      top_sample_clones_change_sample1_df[i, ] = unlist(out[[1]])
      top_sample_clones_change_sample2_df[i, ] = unlist(out[[2]])
      clonal_replacement_df[i, ] =  unlist(out[[3]])
    }
    
  }
  
  outfile = paste0(outdir,"/",compared_pair[[1]], "_vs_", compared_pair[[2]], ".txt" )
  outfile =gsub(" ", "_", outfile)
  write.table(clonal_replacement_df, outfile, sep="\t", quote= FALSE, row.names = TRUE)
  
  return(list(top_sample_clones_change_sample1_df= top_sample_clones_change_sample1_df,
              top_sample_clones_change_sample2_df= top_sample_clones_change_sample2_df,
              clonal_replacement_df= clonal_replacement_df))
}

run_example = function()
{
  
  setwd(
    "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/"
  )
  trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt"
  
  metadata = get_metadata_sorted()
  trackClonotype_df <-
    read.delim(trackClonotype_df_file, check.names = FALSE)
out = get_top_sample_clones_change_clonal_timepoints( compared_pair = c("Blood_Pre-TTF_sc", "Blood_Post-TTF_sc"),
                                                                 metadata= metadata,
                                                                 tops  = c(10, 20, 50, 100),
                                                                 trackClonotype_df = trackClonotype_df)

}