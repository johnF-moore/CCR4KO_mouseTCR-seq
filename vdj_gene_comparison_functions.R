library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)

##-------
# VDJ gene comparison functions for CCR4 mouse TCR-seq project
##-------
get_gene_usage= function(dir_path, md_path, gene= "V"){
  ## Takes iRepertoire raw sequencing data and produces combined data frame of gene usage for all samples
  dirs_out= list.dirs(dir_path, full.names= FALSE, recursive= FALSE)
  
  out_list0= vector('list', length= length(dirs_out))
  out_list1 = out_list0
  
  for(i in 1:length(dirs_out)){
    current_dir= paste0(dir_path,"/", dirs_out[i])
    gene_files= list.files(current_dir,  
                           pattern= paste0(gene,"_\\d_usage.csv"))
    
    gene0= read_csv(paste(current_dir,"/", gene_files[1], sep = ""),
                    col_names= FALSE, show_col_types= FALSE) %>%
      mutate(sample= gsub("_\\w+.csv","", gene_files[1]), 
             origin= as.factor(0)) %>%
      setNames(nm= c("gene", "score", "samples", "origin"))
    
    gene1= read_csv(paste(current_dir,"/",gene_files[2], sep=""),
                    col_names= FALSE, show_col_types= FALSE) %>% 
      mutate(sample= gsub("_\\w+.csv","", gene_files[2]), 
             origin= as.factor(1)) %>%
      setNames(nm= c("gene", "score", "samples", "origin"))
    
    out_list0[[i]]<- gene0
    out_list1[[i]]<- gene1
  }
  cdf <- c(out_list0, out_list1)
  combined_df <- do.call(rbind, cdf)
  combined_df$samples <- paste("s", combined_df$samples, sep= "_")
  
  ## Adding metadata
  metadata <- read_csv(md_path, show_col_types= FALSE)
  
  combined_df <- left_join(combined_df, metadata, by= "samples") %>%
    mutate_all(~replace(., is.na(.), 0))
  combined_df <- combined_df %>% mutate(gene_family = gsub("-\\w+", "", gene))
  
  return(combined_df)
}

score_gene_usage <- function(df, new_name, group_cols="gene", score_col= "score"){
  ## Flexible function to score gene usage across groups
  new_name= sym(new_name)
  df <- df %>% 
    group_by(across(all_of(group_cols))) %>% 
    mutate(sum_score= sum(.data[[score_col]])) %>% 
    rename("{{new_name}}" := sum_score) %>% 
    ungroup()
  
  return(df)
}

# load_genes <- function(path){
#   df <- read_csv(path, show_col_types= FALSE) %>% 
#     mutate(score= as.numeric(score))
#   return(df)
# }

fisher_data_prep <- function(df, gene_col, names_from, values_from){
  fisher_data <- df %>% 
    select(all_of(c(gene_col, names_from, values_from))) %>% 
    distinct() %>% 
    pivot_wider(names_from  = all_of(names_from),
                values_from = all_of(values_from)) %>% 
    column_to_rownames(gene_col) %>% 
    mutate(across(everything(), ~round(.x*10)))
  return(fisher_data)
}

lfc <- function(df, split_col, constant_col, score_col,
                constant_val, pseudo_zero= 1e-6){
  ## Log fold change function for working with formatted differential gene expression data
  uni_split= paste0("log_score_", unique(df[,split_col, drop= TRUE]))
  
  if(length(uni_split) != 2){
    stop("trying to split on a column with more than 2 unique variables.")
  }
  
  lfc <- df %>%
    as.data.frame() %>% 
    select(gene, all_of(c(split_col, constant_col, score_col))) %>% 
    distinct() %>% 
    filter(get(constant_col) == constant_val) %>% 
    ## get() gets columns that match its string argument
    mutate(log_scores= log2(get(score_col) + pseudo_zero)) %>%
    ## have to add a small constant to prevent zeros from going to infinity
    pivot_wider(id_cols= all_of(c("gene", constant_col)),
                values_from= log_scores, 
                names_from= split_col, 
                names_prefix= "log_score_") %>%
    mutate(lfc=  get(uni_split[1]) - get(uni_split[2])) %>% 
    as.data.frame()
  ## Making two columns for the different log values.
  ## This changes the shape of the df.
  ## If you want to add the logFC back to your main df, do pivot_longer and add back with a join
  return(lfc)
}

VDJ_usage_DE_formatter <- function(df, score_col){
  wide_df <- df %>% 
    pivot_wider(id_cols= "gene",
                names_from= c("cell_type","phenotype","samples"), 
                values_from= score_col) %>% 
    rename(X= gene)
  
  colnames(wide_df) <- gsub("s_","", colnames(wide_df))
  
  return(wide_df)
}

widyr_wilcox <- function(df, group1, group2){
  rows= nrow(df)
  p_vals= vector(mode= "numeric",length= rows)
  
  for(row in 1:rows){
    x= df[row,] %>% select(contains(group1)) %>% as.numeric()
    y= df[row,] %>% select(contains(group2)) %>% as.numeric()
    p_vals[row]= suppressWarnings(wilcox.test(x= x, y= y)$p.value)
  }
  p_vals= p.adjust(p= p_vals, method= "BH", n= rows)
  
  return(p_vals)
}

widyr_fisher <- function(df, group1, group2){
  rows= nrow(df)
  p_vals= vector(mode= "numeric", length= rows)
  for(row in 1:rows){
    x= df[row,] %>% select(contains(group1)) %>% as.numeric()
    y= df[row,] %>% select(contains(group2)) %>% as.numeric()
    if(length(unique(x)) == 1 || length(unique(y)) == 1){
      p_vals[row]= NA
      ## fisher.test doesn't work if there aren't 2 unique values for each group
    } else{
      p_vals[row]= fisher.test(x= x, y= y, simulate.p.value = TRUE, B= 2000)$p.value
    }
  }
  p_vals= p.adjust(p= p_vals, method= "BH", n= rows)
  return(p_vals)
}

get_gene_reads= function(dir_path, md_path){
  ## Gets all of the raw data that we need for differential expression from all TCR-seq samples 
  dirs_out= list.dirs(dir_path, full.names= FALSE, recursive= FALSE)
  
  out_list= vector('list', length= length(dirs_out))
  
  for(i in 1:length(dirs_out)){
    current_dir= paste0(dir_path,"/", dirs_out[i])
    gene_file= list.files(current_dir, pattern= "_CDR3_list_2.csv")
    
    gene_reads= read_csv(paste(current_dir,"/", gene_file, sep = ""),
                         col_names= FALSE, show_col_types= FALSE) %>%
      setNames(nm= c("X", "v_gene", "j_gene", "reads")) %>% 
      mutate(samples= gsub("_\\w+.csv","", gene_file))
    
    out_list[[i]]<- gene_reads
  }
  combined_df <- do.call(rbind, out_list)
  combined_df$samples <- paste("s", combined_df$samples, sep= "_")
  
  ## Adding metadata
  metadata <- read_csv(md_path, show_col_types= FALSE)
  
  combined_df <- left_join(combined_df, metadata, by= "samples") %>%
    mutate_all(~replace(., is.na(.), 0))
  combined_df <- combined_df %>% mutate(v_subgroup = gsub("(\\w+\\d+).*","\\1", v_gene))
  ## There are no IMGT mTRAJ subgroups.
  return(combined_df)
}

