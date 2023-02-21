###Install the required libraries-------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
#--------------------------------------------


##Clear environment------------------
rm(list=ls()) 
#------------------------------------



###Load the required libraries-------------
library("DESeq2")
library("edgeR")
#------------------------------------------




###constants-------------------------------
directory <- ("CodeDirectory")

#the names of the data sets [will be used to load the count and meta data files]
#datasets should be named as [ds_name]_data & [ds_name]_meta in the directory above
#exmple: CKD-AmericanHuman_data.csv & CKD-AmericanHuman_meta.csv
data_sets <- c("CKD-AmericanHuman", "CKD-ChineseHuman", "CKD-FP-mouse") #we created another mouse ds to relevel

#control groups for the different data sets--if we set intercept <- T
control_group <- c("Control", "Control", "Sham", "CKD")

#set the contrasts that we need:
contrasts <- list(
  list(list("Group", "CKD", "Control")),
  list(list("Group", "CKD", "Control"), list("Group", "CKD_HTN", "Control")),
  list(list("Group", "Sham_FP", "Sham"), list("Group", "CKD", "Sham"), list("Group", "CKD_FP", "CKD"))
)


#set the working directory
setwd(directory)

#define if we need reference group or not-will affect creating the contrasts
#if TRUE the groups in the control_group list will be used as reference for the data sets
#We will not use intercept, as we will need two references for the mouse DS
intercept <- F
#-------------------------------------------


###Here we iterate over the datasets, load the count matrices and meta data then do the DE------------------------------------
for (ds in 1:length(data_sets))
{
  message(paste("performing DGE analysis for dataset:", data_sets[ds]))
  
  ###reading files----------------------------
  message("creating the name of the files...")
  count_file <- paste0(data_sets[ds], "_data.csv")
  meta_file <- paste0(data_sets[ds], "_meta.csv")
  
  message(paste("loading the count matrix file:", count_file))
  count_data <- as.data.frame(read.csv(count_file, row.names = "feature", check.names = F))
 
  message(paste("loading the metadata file:", meta_file))
  col_data <- as.data.frame(read.csv(meta_file, row.names = "sample", check.names = F))
  col_data$Group <- as.factor(col_data$Group)

  #check if the column names are consistent [must be TRUE]
  if(!all(rownames(col_data) == colnames(count_data)))
  {
    stop("Count matrix's columns' names and meta data columns' names are not consistent!")
  }
  #-------------------------------------------
  
  
  ###filter readings---------------------------
  message("clean up the data...")
  #create a temporary DGEList object to filter the readings
  temp_dge<- DGEList(count_data)
  temp_dge$samples$group <- as.factor(col_data$Group)
  goodExpressions <- filterByExpr(temp_dge, group='group')
  temp_dge <- temp_dge[goodExpressions,, keep.lib.sizes=FALSE]
  #save the filtered count to the count_data object
  count_data <- temp_dge$counts
  #free the memory from the DGEList object
  rm(temp_dge)
  write.csv(count_data, file=paste0("./", data_sets[ds], "/filtered_counts.csv"))
  #-------------------------------------------
  
  
  
  ###perform differential analysis------------------
  #NOTE:DESeq2 doesn't use normalized counts, rather, it uses the raw counts and models the normalization 
  #     within the Generalized Linear Model (GLM)-negative binomial distribution 
  message("creating the DESeqDataSet object...")
  
  
  #create the design matrix based on the value of reference variable [with ref group or not]
  if (intercept)
  {
    dds <- DESeqDataSetFromMatrix (
      countData = count_data, 
      colData = col_data,
      design = ~ Group
    )
    #set the reference level
    dds$Group <- relevel(dds$Group, ref = control_group[ds])
  }
  else
  {
    dds <- DESeqDataSetFromMatrix (
      countData = count_data, 
      colData = col_data,
      design = ~-1 +  Group
    )
  }
  
  
  message("fitting the GLM model...")
  dds <- DESeq(dds, sfType = 'poscounts')
  
  #create a folder for each data set
  if(!dir.exists(data_sets[ds]))
    dir.create(data_sets[ds])
  

  #get the normalized counts
  normalized_counts <- counts(dds, normalized=TRUE)
  write.csv(normalized_counts, file=paste0("./", data_sets[ds], "/normalized counts.csv"))
  

  
  for (ds_cont in 1:length(contrasts[[ds]]))
  {
    current_contrast <- unlist(contrasts[[ds]][[ds_cont]])
    #friendly format for the contrast
    contrast_name <- paste0(current_contrast[2], "_vs_", current_contrast[3]) 
   
    #get the diff table for the contrast
    res <- results(dds, contrast = current_contrast, alpha = 0.05)
    message(paste("summary results for:", contrast_name))
    summary(res)
    
    #remove NA padj rows
    res <- res[!is.na(res$padj), ]
    
    #order the results based on the adjusted p-value
    res <- res[order(res$padj), ] 
    
    #combine the normalized counts with the DEG statistics-for convenience
    out <- cbind(as.data.frame(res), normalized_counts[rownames(res), ])
    write.csv(out, file=paste0("./", data_sets[ds], "/", contrast_name, "_de_ALL.csv"))
  }
  #-------------------------------------------------
}
###-----------------------------------------------------------------------------------------------------------------------------
