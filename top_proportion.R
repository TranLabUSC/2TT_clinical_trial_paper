trackClonotype_df <- read.delim("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt", check.names = FALSE)
library(dplyr)
get_top_proportion= function(top = 10, trackClonotype_df)
{
  df = trackClonotype_df[,2:ncol(trackClonotype_df)]
  top_proportion = sapply(df,function(x)
    {
    x_sorted = sort(x, decreasing = TRUE)
    top_prop = sum(x_sorted[1:top])
    return(top_prop)
  })
  top_proportion = as.data.frame(top_proportion)
  colnames(top_proportion) = paste0("top_", top, "_proportion")
  return(top_proportion)
}

make_top_proportion_table = function(trackClonotype_df_file = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt" ,
                                     tops = c(10,20,50,100),
                                     outfile= "top_proportions.txt")
{
  
  trackClonotype_df <- read.delim(trackClonotype_df_file, check.names = FALSE)
  top_prop_table = NULL
  for (top in tops)
  {
    top_prop = get_top_proportion(top = top, trackClonotype_df = trackClonotype_df)
    if (is.null(top_prop_table))
    {
      top_prop_table = top_prop
    }else
    {
      top_prop_table = cbind(top_prop_table, top_prop)
    }
  }
  write.table(top_prop_table, file= outfile, sep="\t", quote = FALSE, row.names=TRUE)
  return(top_prop_table)
}

make_top_proportion_table()

