#install.packages("dplyr")
library(dplyr)


get_cell_fraction_from_scRNAseq = function(filtered_contig_annotations_file = "filtered_contig_annotations.csv")
{
  setwd(
    "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/2TT_aggregate/outs/vdj_t"
  )
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

# Now visualization

library(immunarch)



#dev.off()
#metadata_preparation
prepare_metadata = function()
{
  # Single cell meta
  sc_meta <- read.csv("Metadata/scRNAseq_meta.csv")
  sc_meta$timepoint = lapply(sc_meta$origin, function(x) {
    strsplit(x, "_")[[1]][1]
  })
  sc_meta[sc_meta$timepoint == "N", "timepoint"] = "Pre-TTF"
  sc_meta[sc_meta$timepoint == "Y", "timepoint"] = "Post-TTF"
  sc_meta[sc_meta$timepoint == "T", "timepoint"] = "Tumor"
  sc_meta$timepoint = factor(
    sc_meta$timepoint,
    levels = c (
      "Tumor",
      "Pre-TTF",
      "Post-TTF",
      "1",
      "2",
      "3",
      "4",
      "5",
      "R",
      "R2"
    )
  )
  sc_meta = sc_meta[with(sc_meta, order(Patient, timepoint)), ]
  sc_meta = sc_meta[, c("sample_id", "Patient", "Sample.Type", "timepoint")]
  sc_meta$dataset = "sc"
  sc_meta$sample_id = paste0(sc_meta$sample_id, "_sc")
  sc_meta$Patient = paste0("p", sc_meta$Patient)
  
  
  # Immunoseq meta
  immunoseq_meta <- read.delim("Metadata/immunoseq_meta.csv")
  immunoseq_meta = immunoseq_meta %>% filter(Trial == '2TT')
  immunoseq_meta = immunoseq_meta[, c("Sample.Name", "Patient_ID", "Sample.type", "TimePoint")]
  colnames(immunoseq_meta) = c("sample_id", "Patient", "Sample.Type", "timepoint")
  immunoseq_meta$Patient = paste0("p", immunoseq_meta$Patient)
  immunoseq_meta$dataset = "immunoseq"
  
  #bulk_meta
  bulk_meta <-
    read.delim(
      "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/Metadata/bulkRNAseq_meta.csv"
    )
  bulk_meta$sample_type = "Blood"
  bulk_meta[bulk_meta$group == "T_0", "timepoint"] = "Pre_TTF"
  bulk_meta[bulk_meta$group == "T_1", "timepoint"] = "Post_TTF"
  bulk_meta = bulk_meta[, c("sample", "patient", "sample_type", "timepoint")]
  
  colnames(bulk_meta) =  c("sample_id", "Patient", "Sample.Type", "timepoint")
  bulk_meta$dataset = "bulk"
  bulk_meta$sample_id = paste0("bulk_", bulk_meta$sample_id)
  combined_meta = rbind(sc_meta, immunoseq_meta, bulk_meta)
  write.table(
    combined_meta,
    file = "Metadata/combined_meta.txt" ,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  return(combined_meta)
}




track_patient_clone = function(patient, top, filename, metadata)
{
  out = tryCatch({
    patient = "p18"
    patient_samples = metadata %>% filter(Patient == patient)
    
    #tc1  =trackClonotypes(clonotyping_data$data[patient_samples$sample_id], list(nrow(patient_samples), 1000), .col = "aa")
    
    tc1  = trackClonotypes(clonotyping_data$data[patient_samples$sample_id], list("18P63_sc", 10000), .col = "aa")
    #tc2  =trackClonotypes(clonotyping_data$data[patient_samples$sample_id], list("T24F_TCRB", 10000), .col = "aa")
    #tc3  =trackClonotypes(clonotyping_data$data[patient_samples$sample_id], list("bulk_29_24-1", 10000), .col = "aa")
    
    tc1  = trackClonotypes(clonotyping_data$data[c(patient_samples$sample_id,
                                                   "bulk_29_24-1",
                                                   "bulk_30_24-2")], list(1, 10000), .col = "aa")
    tc2  = trackClonotypes(clonotyping_data$data[patient_samples$sample_id], list("25P63_sc", 10000), .col = "aa")
    tc3  = trackClonotypes(clonotyping_data$data[patient_samples$sample_id],
                           list("bulk_31_25-1", 10000),
                           .col = "aa")
    dt = clonotyping_data$data[["bulk_31_25-1"]]
    
    
    
    
    
    
    png(
      file = filename ,
      units = "px",
      width = 250 * nrow(patient_samples),
      height = 1000
    )
    p = vis(tc1) + theme_bw() +
      theme(
        legend.position = "none" ,
        text = element_text(size = 30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 40)
      ) +  scale_y_continuous(limits = c(0, 1))  +
      scale_x_discrete(labels =  patient_samples$timepoint)  +
      annotate(
        geom = "text",
        x = 1,
        y = 0.9,
        label = paste("P", patient),
        size = 20
      )
    print(p)
    dev.off()
  }
  , error = function(cond)
  {
    print(cond)
  })
  
  
}

# track all clones by aa :

get_all_unique_sequences = function(.col = "aa",
                                    clonotyping_data = clonotyping_data)
{
  # .col can be "aa" or "nt": amino acid or nucleotide
  select_col = paste0("CDR3.", .col)
  sequences = NULL
  for (df in clonotyping_data$data)
  {
    df = as.data.frame(df)
    if (!"chain" %in% colnames(df))
    {
      df$chain = "TRB"
    }
    sample_sequences = df[, c(select_col, "chain")]
    if (is.null(sequences))
    {
      sequences = sample_sequences
    } else
    {
      sequences = rbind(sequences, sample_sequences)
    }
  }
  sequences = sequences[!duplicated(sequences), ]
  return(sequences)
}


#combine all aa CDR3 sequence
#tc_all  = trackClonotypes(clonotyping_data$data, .col = "aa")


get_TRB_only_data = function(df)
  #This function correct single cell clonotyping data, removing TRA sequence . If the row does not contain TRB sequence, remove the row. At the end, the fraction is recalculated.
  #If the row contain two TRB chain, then split into two clones with equal proportion. 
{
  new_df = NULL
  for (i in c(1:nrow(df)))
  {
    row = df[i, ]
    chains = as.character(row["chain"])
    chains = unlist(strsplit(chains, split = ";"))
    for (j in 1:length(chains))
    {
      if (chains[j] == "TRB")
      {
        if (is.null(new_df))
        {
          new_df = row
        } else
        {
          new_df = rbind(new_df, row)
        }
        
      }
    }
    
  }
  
}


trackClonotypesRelaxedAll = function(clonotyping_data, .col = "aa")
  # Cells belong to the same clone if ANY of the chain sequences are the same, not ALL as in ImmunArch. This is relevant to scRNAseq as some single cells only contain sequences only for TRA, or TRB, while other contain TRA, TRA, TRB, TRB due to leakiness of allelic exclusion.
  #.col. can be aa or nt. (Amino acid or nucleotide)
  # Algorithm:
  #1)collect all unique sequence combinations
  #2) For each unique sequence combinations:
  # Check if the clone already exist:
  # If exist: update the clone, sequence combination.
  # If not: create new clone
  # Clone structure:
  #cloneID
  #TRA sequences
#TRB sequences
#combination variations
# Result: not good, too many cells with many combinations. Those most likely doublet.
{
  #1)collect all unique sequence combinations
  
  sequence_combinations = get_all_unique_sequences(.col = .col, clonotyping_data = clonotyping_data)
  
  
  
  #2) For each unique sequence combinations:
  # Check if the clone already exist: using named vector(key is chain_sequence)
  # If exist: update the clone, sequence combination.
  # If not: create new clone
  # Clone structure:
  #cloneID
  #TRA sequences: vector
  #TRB sequences: vector
  #combination variations: dataframe columns: sequence, chain.
  
  # Create clones from sequence combination
  
  clones = get_clone_from_sequence_combination(sequence_combinations)
  
}


get_clone_from_sequence_combination = function(sequence_combinations)
  #sequence combination is a data frame with two columns first one is sequence, second is chains
{
  clones = list()
  clone_id = 1
  clone_dict = vector()
  for (i in 1:nrow(sequence_combinations))
  #for (i in 1:80000)
  {
    combination = sequence_combinations[i, ]
    chains = as.character(combination["chain"])
    chains = unlist(strsplit(chains, split = ";"))
    combination_sequences = as.character(combination[1])
    combination_sequences =  unlist(strsplit(combination_sequences, split = ";"))
    #Create first clone
    if (length(clones) == 0)
    {
      out = create_new_clone(clone_id, clone_dict, chains, combination_sequences,combination)
      clone = out$clone
      clones[[clone_id]] = clone
      
      clone_dict = out$clone_dict
      clone_id = clone_id + 1
     
    } else
      
    {
      clone_exists = FALSE
      for (j in 1:length(chains))
      {
        chain = chains[j]
        sequence_ = combination_sequences[j]
        key = paste0(chain, "_", sequence_)
        # Clone exists but different combination, key exists
        if (!is.na(clone_dict[key]))
        {
          clone_exists = TRUE
          existing_clone_id = as.integer(clone_dict[[key]])
          #if (existing_clone_id  == 64560)
          #{
          #  print(paste("clone exit, updating clone. Clone ID is:", clone_dict[key], "Clone key is:", key, "iteration:", i))
          #}
          clone = clones[[existing_clone_id]]
          clone[["combinations"]] = rbind(clone[["combinations"]] , as.data.frame(combination))
          #print("clone combinations:")
          #print(clone[["combinations"]])
          
          # update TRA, TRB
          for (k in 1:length(chains))
          {
            sequence_ = combination_sequences[k]
            chain = chains[k]
            if (chain == "TRA" |  chain == "TRB" )
            {
              if (!sequence_ %in% clone[[chain]])
              {
                clone[[chain]] = c(clone[[chain]], sequence_)
                clone_dict[paste0(chain, "_", sequence_)] = existing_clone_id
               
              }
              
            } else
            {
              stop("The chain is not TRA, TRB, only TRA or TRB or accepted")
            }
            
          }
          
          # Replace the old clone with the updated one
          clones[[existing_clone_id]] = clone 
          break
        }
      }
      # if no of the sequence exists, then create new clone
      if (!clone_exists)
      {
        out = create_new_clone(clone_id, clone_dict, chains, combination_sequences,combination)
        clone = out$clone
        clones[[clone_id]] = clone
        clone_dict = out$clone_dict
        clone_id = clone_id + 1
       
      }
    }
  }
  
  n_comb_df_all = NULL
  for (clone in clones)
  {
    clone_id = clone[["clone_id"]]
    print(clone_id)
    n_comb = nrow(clone[["combinations"]])
    n_comb_df = data.frame(clone_id = clone_id, n_comb = n_comb)
    if (is.null(n_comb_df_all))
    {
      n_comb_df_all = n_comb
    }else
    {
      n_comb_df_all = rbind(n_comb_df_all, n_comb_df)
    }
  }
  
  n_comb_df_all = n_comb_df_all[order(n_comb_df_all$n_comb, decreasing = TRUE), ]
  
  return(list(clones, clone_dict, n_comb_df_all))
}


create_new_clone = function(clone_id, clone_dict,
                            chains,
                            combination_sequences, combination)
{
  # get TRA, TRB
  TRA = vector()
  TRB = vector()
  clone_id = as.integer(clone_id) # Make sure that clone_id is integer, not a data frame extract
  for (j in 1:length(chains))
  {
    sequence_ = combination_sequences[j]
    chain = chains[j]
    if (chain == "TRA")
    {
      TRA = c(TRA, sequence_)
    } else if (chain == "TRB")
    {
      TRB = c(TRB, sequence_)
    } else
    {
      stop("The chain is not TRA, TRB, only TRA or TRB or accepted")
    }
    clone_dict[paste0(chain, "_", sequence_)] = clone_id
  }
  
  clone = list(
    clone_id = clone_id,
    TRA = TRA,
    TRB = TRB,
    combinations = as.data.frame(combination)
  )
  print(paste("Made new clone  clone id made:", clone_id))
  return(list(clone = clone, clone_dict = clone_dict))
}

get_frequency = function(df)
{
  chain_proportion =  as.data.frame(table(df[, "chain"]))
  chain_proportion$relative = chain_proportion$Freq / sum(chain_proportion$Freq)
  write.table(
    chain_proportion,
    file = "scT_chain_proportion.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  return(chain_proportion)
}


get_result = function(tops = c(10, 50, 100, 200))
{
  image_folder = "clone_tracking_endpoint"
  dir.create(image_folder)
  for (top in tops)
  {
    dir.create(paste0(image_folder, "/top", top))
    
  }
  
  
  for (patient in patients)
  {
    for (top in tops)
    {
      print(paste(patient, "top", top))
      filename = paste0(image_folder, "/top", top, "/P", patient, ".png")
      track_patient_clone(patient, top, filename)
    }
    
    
  }
  
}

#setwd(
#  "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping"
#)

#data_folder = "~/Desktop/USC_Cluster/Data/Collaborator/DJ/scRNAseq/latest_big_sequencing_projects/2TT/Analysis/Aggregate/Clonotyping/ClonotypingData"

#clonotyping_data = repLoad(data_folder)

#load(".RData")
#.col = "aa"
#sequence_combinations = get_all_unique_sequences(.col = .col, clonotyping_data = clonotyping_data)
#out =  get_clone_from_sequence_combination(sequence_combinations)
#trackClonotypesRelaxedAll( clonotyping_data = clonotyping_data)

#clonotyping_data_sc  = repLoad(data_folder)
