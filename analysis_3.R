setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/code"
)
source("analysis.R")

setwd(
  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
)

data_folder = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/ClonotypingData"

#clonotyping_data = repLoad(data_folder)

#load(".RData")
.col = "aa"
#sequence_combinations = get_all_unique_sequences(.col = .col, clonotyping_data = clonotyping_data)

# get_clone_dict = function(sequence_combinations)
# {
#   clone_id = 1
#   clone_dict = vector()
#   duplicated_clone = vector()
#   for (i in 1:nrow(sequence_combinations))
#     # only consider TRB . One TRB is one clone. If in the data, clone contain more than 1 TRB, then evenly split int subclones where each sublcone has 1 TRB.
#   {
#     print(i)
#     combination = sequence_combinations[i,]
#     chains = as.character(combination["chain"])
#     chains = unlist(strsplit(chains, split = ";"))
#     combination_sequences = as.character(combination[1])
#     combination_sequences =  unlist(strsplit(combination_sequences, split = ";"))
#     duplicated_clone = vector
#     #Create first clone
#     for (j in 1:length(chains))
#     {
#       chain = chains[j]
#       if (chain == "TRB")
#       {
#         sequence_ = combination_sequences[j]
#         if (is.na(clone_dict[sequence_]))
#         {
#           clone_dict[sequence_] = clone_id
#           clone_id = clone_id + 1
#           
#         } else
#         {
#           print(paste("sequence is already in the dictionary", combination))
#           duplicated_clone = c(duplicated_clone, combination)
#         }
#       }
#       
#     }
#     
#   }
#   
#   return(clone_dict)
# }
# 
# 
# get_unique_TRB_sequences = function(sequence_combinations)
# {
#   TRB_sequences = vector()
#   for (i in 1:nrow(sequence_combinations))
#     # only consider TRB . One TRB is one clone. If in the data, clone contain more than 1 TRB, then evenly split int subclones where each sublcone has 1 TRB.
#   {
#     print(i)
#     combination = sequence_combinations[i,]
#     chains = as.character(combination["chain"])
#     chains = unlist(strsplit(chains, split = ";"))
#     combination_sequences = as.character(combination[1])
#     combination_sequences =  unlist(strsplit(combination_sequences, split = ";"))
#     
#     #Create first clone
#     for (j in 1:length(chains))
#     {
#       chain = chains[j]
#       if (chain == "TRB")
#       {
#         sequence_ = combination_sequences[j]
#         
#         TRB_sequences = c(TRB_sequences, sequence_)
#       }
#     }
#     
#     
#   }
#   
#   
# }


trackClonotypes_TRB = function(clonotyping_data, .col = "aa")
  
{
  select_col = paste0("CDR3.", .col)
  trackClonotypes_df = data.frame()
  for (n in 1:length(clonotyping_data$data))
  #for (n in 1:3)
  {
    sample = names(clonotyping_data$data)[[n]]
    print(paste(n, sample))
    
    sample_data = clonotyping_data$data[[sample]]
    # parsing the chain
    # Loop through each clone in the data
    for (i in 1:nrow(sample_data))
    {
      n_clones = sample_data[[i, "Clones"]]
      
      V = sample_data[[i, "V.name"]]
      sequences = sample_data[[i, select_col]]
      V = unlist(strsplit(V, split = ";"))
      sequences = unlist(strsplit(sequences, split = ";"))
      # loop through each chain
      for (j in 1:length(V))
      {
        if (substr(V[j], 1, 3) == "TRB")
        {
          if (is.null(trackClonotypes_df[sequences[[j]], sample])  )
              {
            trackClonotypes_df[sequences[[j]], sample] = n_clones
          }else if ( is.na(trackClonotypes_df[sequences[[j]], sample]) )
          {
            trackClonotypes_df[sequences[[j]], sample] = n_clones
            
          }else
                {
                  trackClonotypes_df[sequences[[j]], sample] = trackClonotypes_df[sequences[[j]], sample] +  n_clones
                }
          
        }
      }
    }
    
    
  }
  
  trackClonotypes_df[,select_col] = row.names(trackClonotypes_df)
  trackClonotypes_df =  trackClonotypes_df[, c(ncol(trackClonotypes_df), 1:(ncol(trackClonotypes_df)-1) )]
  write.table(trackClonotypes_df, file="trackClonotype_df.abs.txt", sep="\t", quote=FALSE)
  trackClonotypes_df[is.na(trackClonotypes_df)] = 0
  
  trackClonotypes_df_relative = trackClonotypes_df
  trackClonotypes_df_relative[,2:ncol(trackClonotypes_df)] = sapply(trackClonotypes_df_relative[,2:ncol(trackClonotypes_df)] , function(x)
  {
    x/sum(x)
  })
  
  write.table(trackClonotypes_df_relative, file="trackClonotype_df.proportion.txt", sep="\t", quote=FALSE)
  return(trackClonotypes_df)
}

trackClonotypes_df = trackClonotypes_TRB(clonotyping_data = clonotyping_data)

trackClonotypes_df




#tc1  = trackClonotypes(clonotyping_data$data[1:4], list(1, 200), .col = "aa")

# Clone align directly

