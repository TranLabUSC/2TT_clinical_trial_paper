

# for each vector:
 simpson_diversity_index = function(clonal_proportions)
 {
   
   # Calculate
   #Simpson's Diversity Index (SDI) 
   #SDI = 1 - ∑(ni / N)^2
   # a vector, where each value is the proportion of each clone. Total proportion is 1
   # Gini-Simpson
   SDI = 1 - sum(clonal_proportions^2)
    return(SDI)  
 }
 
 shannon_diversity_index = function(clonal_proportions)
 {
   
   # Calculate
  
   #H = - sum(p*ln(pi))
   clonal_proportions = clonal_proportions[clonal_proportions > 0]
   h = (clonal_proportions * log(clonal_proportions))
   h = -sum(h)
   return(h)  
 }
 

get_clonotype_diversity = function(trackClonotype_df.relative, diversity_index = "shannon", outfolder=".")
{
  diversity_dict = list(shannon = shannon_diversity_index,   simpson = simpson_diversity_index )
  
clonotype_diversity = sapply(trackClonotypes_df_relative[,2:ncol(trackClonotypes_df_relative)],diversity_dict[[diversity_index]]
  )
clonotype_diversity = as.data.frame(clonotype_diversity)
write.table(clonotype_diversity, paste0(outfolder, "/", diversity_index, "_clonal_diversity.txt"), sep="\t", row.names=TRUE, quote=FALSE)
return(clonotype_diversity)
}

run = function()
# setwd("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping")
{trackClonotypes_df_relative <- read.delim("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt", check.names = FALSE)
shannon = diversity = get_clonotype_diversity(trackClonotype_df.relative = trackClonotype_df_relative, diversity_index = "shannon")
simpson =diversity = get_clonotype_diversity(trackClonotype_df.relative = trackClonotype_df_relative, diversity_index = "simpson")
return(list(shannon = shannon, simpson = simpson))
}

 out = run()
 print(out[[1]])
 print(out[[2]])