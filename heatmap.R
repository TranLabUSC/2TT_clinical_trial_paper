#heat map analysis
setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)
trackClonotype_df <- read.delim("~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/trackClonotype_df.proportion.txt")
library(corrplot)
corr = cor(trackClonotype_df[,2:ncol(trackClonotype_df)])
corrplot(corr,order="hclust", hclust.method = "average", addgrid.col = NA, diag=FALSE, tl.cex = 0.12, type= "upper", tl.col="black", addrect = 20,  col=colorRampPalette(c("blue","white","red"))(200))
