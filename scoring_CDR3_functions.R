## John Moore
## 08/31/22
## Functions for scoring CDR3 using PWMs
library(ggplot2)
library(dplyr)
library(tibble)
library(DT)
library(here)
library(tidyr)
library(data.table)
library(ggpubr)
library(reactable)

##--------
# Functions
##--------
PWM_list_to_df= function(PWMs){ ## Where PWMs is a nested list 
  for(i in 1:length(PWMs)){
    pheno.list= PWMs[[i]]
    CDR3_lens= as.numeric(names(pheno.list))
    
    for(j in 1:length(pheno.list)){
      pheno.list[[j]]$CDR3_len= CDR3_lens[j]
      pheno.list[[j]]= pheno.list[[j]] %>% 
        pivot_wider(names_from= aa, values_from= prop) %>% 
        as.data.frame()
    }
    
    PWMs[[i]]= do.call(rbind, pheno.list)
    PWMs[[i]]$phenotype= names(PWMs)[i]
    
  }
  PWMs= do.call(rbind, PWMs)
  return(PWMs)
}


load_tidy_PWM= function(path){
  PWMs= readRDS(path) %>% 
    PWM_list_to_df() %>% 
    ## Convert nested lists to df
    setDT() %>%
    melt(id.vars= c("position","CDR3_len","phenotype"), 
         variable.name= "aa", 
         value.name= "prop") %>% 
    as.data.frame()
  ## Melt to tidy 
  
  return(PWMs)
}


convert_to_aa_group= function(df, col){
  col_expr= expr(col)
  side_chains= df %>% 
    mutate(
      {{col_expr}} := case_when(
        grepl("A|P|G|I|L|V", .data[[col]]) ~ "Alipathic",
        grepl("F|W|Y", .data[[col]]) ~ "Aromatic",
        grepl("D|E", .data[[col]]) ~ "Acidic",
        grepl("R|H|K", .data[[col]]) ~ "Basic",
        grepl("S|T", .data[[col]]) ~ "Hydroxylic",
        grepl("N|Q", .data[[col]]) ~ "Amidic",
        grepl("C|M", .data[[col]]) ~ "Sulfur-containing",
        grepl("*", .data[[col]]) ~ "*"
      ) 
    ) %>% 
    pull({{col_expr}})
  
  return(side_chains)
}


score_PWM= function(ex_CDR3, PWM, is_side_chain, replace=0.0001){
  ## Checking inputs
  if(suppressWarnings(!is.na(as.numeric(ex_CDR3)))){
    stop("CDR3 is not a character or factor.")
  } else if(is.factor(ex_CDR3)){
    ex_CDR3= as.character(ex_CDR3)
  }
  if(!is.data.frame(PWM) && is.matrix(PWM)){
    PWM= as.data.frame(PWM)
  }

  ## Splitting amino acids & setting up filtering
  aa_CDR3= strsplit(ex_CDR3,"")[[1]]
  
  if(is_side_chain){
    aa_CDR3= aa_CDR3 %>% 
      as.data.frame() %>% 
      convert_to_aa_group(col= names(.))
  }
  aa_pos_CDR3= paste0(aa_CDR3,"_p", 1:length(aa_CDR3))
  indiv_len= nchar(ex_CDR3)
  
  ## Calculating scores for all phenotypes 
  prob_scores= PWM %>% 
    filter(CDR3_len == indiv_len) %>% 
    ## Filtering out rows of incorrect length
    mutate(aa_pos= paste0(aa, "_", position)) %>%
    filter(aa_pos == aa_pos_CDR3) %>% 
    ## Keeping only correct positions and aa
    mutate(prop= ifelse(prop == 0, replace, prop)) %>% 
    ## Replacing 0 w/ some small pseudominimum
    ## I don't want to lose any TCRs just b/c they have a column
    ## that doesn't have that aa.
    ## I'm doing this for our data b/c I don't think our PWM 
    ## is exhaustive. 
    group_by(phenotype) %>% 
    mutate(pheno_score= prod(prop),
           CDR3_aa= ex_CDR3) %>% 
    distinct(phenotype, pheno_score) 
  
  scores= prob_scores$pheno_score
  names(scores)= prob_scores$phenotype
  
  return(scores)
}


score_all_CDR3= function(CDR3, PWM, is_side_chain){
  all_scores= sapply(CDR3, FUN= score_PWM, PWM= PWM, is_side_chain= is_side_chain) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(max= pmax(tconv_WT, tconv_CCR4KO,
                     treg_WT, treg_CCR4KO)) %>% 
    ## Finding max score for each CDR3 
    arrange(desc(max)) %>% 
    rownames_to_column("CDR3_aa")
  return(all_scores)
  ## Output isn't tidy
  ## I may convert to tidy for graphing later, 
  ## but I like this better for table visualization.
}


lm_summary= function(df, x_val, y_val){
  model= lm(df[[y_val]] ~ df[[x_val]])
  plot(df[[y_val]] ~ df[[x_val]], 
       xlab= x_val,
       ylab= y_val)
  abline(model)
  
  return(summary(model))
}





PWM_scoring_wrapper= function(pwm_path, treg_DE_path, tconv_DE_path, is_side_chain){
  ##-----------
  # Loading data
  ##-----------
  PWM= load_tidy_PWM(pwm_path)
  unfilt_treg_TCRs_df= read.csv(treg_DE_path) 
  treg_TCRs_df= unfilt_treg_TCRs_df %>% filter(FDR < 0.99 | abs(logFC) > 2) 
  tconv_TCRs_df= read.csv(tconv_DE_path) %>% filter(FDR < 0.99 | abs(logFC) > 2) 
  
  ##---------
  # Scoring TCRs
  ##---------
  treg_TCR_scores = score_all_CDR3(CDR3= treg_TCRs_df$X,  
                                   PWM= PWM,
                                   is_side_chain = is_side_chain)
  tconv_TCR_scores= score_all_CDR3(CDR3= tconv_TCRs_df$X, 
                                   PWM= PWM,
                                   is_side_chain= is_side_chain) 
  
  ##----------
  # Making lm data
  ##----------
  treg_scores_df= inner_join(treg_TCRs_df,
                             treg_TCR_scores, 
                             by= c("X" = "CDR3_aa")) %>% 
    mutate(treg_max= pmax(treg_WT, treg_CCR4KO),
           abs_logFC= abs(logFC))
  ## Calculating max by cell type for linear modeling
  
  tconv_scores_df= inner_join(tconv_TCRs_df,
                              tconv_TCR_scores, 
                              by= c("X" = "CDR3_aa")) %>% 
    mutate(tconv_max= pmax(tconv_WT, tconv_CCR4KO),
           abs_logFC= abs(logFC))
  
  ##----------
  # CDR3 Tables for referencing DE TCRs 
  ##----------
  treg_DT= treg_scores_df %>% 
    select(-`F`,-PValue, -isRed, -max) %>%
    datatable(caption= "DE and/or High logFC Treg TCRs")
  tconv_DT= tconv_scores_df %>% 
    select(-`F`,-PValue, -isRed, -max) %>%
    datatable(caption= "DE and/or High logFC Tconv TCRs")
  DT_list= list(treg_DT, tconv_DT)
  names(DT_list)= c("treg","tconv")
  
  ##----------
  # PWM validation and modeling
  ##----------
  ## By FDR
  lm_summary(treg_scores_df,  x_val= "FDR", y_val= "treg_max")
  lm_summary(tconv_scores_df, x_val= "FDR", y_val= "tconv_max")
  
  ## By logFC
  lm_summary(treg_scores_df,  x_val= "logFC", y_val= "treg_max")
  lm_summary(tconv_scores_df, x_val= "logFC", y_val= "tconv_max")
  
  ## By abs_logFC
  lm_summary(treg_scores_df,  x_val= "abs_logFC", y_val= "treg_max")
  lm_summary(tconv_scores_df, x_val= "abs_logFC", y_val= "tconv_max")
  
  ##-----------
  # Does probability score differ b/w DE and non-DE TCRs?
  ##-----------
  treg_TCR_scores = score_all_CDR3(CDR3= unfilt_treg_TCRs_df$X, 
                                   PWM= PWM,
                                   is_side_chain= is_side_chain)
  treg_scores_df= inner_join(treg_TCRs_df,
                             treg_TCR_scores, 
                             by= c("X" = "CDR3_aa")) %>% 
    mutate(treg_max= pmax(treg_WT, treg_CCR4KO),
           abs_logFC= abs(logFC))
  
  
  wilcox_pwmScore= compare_means(treg_max ~ isRed, 
                                 data= treg_scores_df, 
                                 method= "wilcox.test",
                                 p.adjust.method= "holm") %>% 
    mutate(y.position=0.0000025)
  
  p= ggplot(treg_scores_df, aes(x= isRed, y= treg_max)) + 
    geom_boxplot() +
    geom_point(position= position_jitter(width= 0.1)) + 
    theme_bw() + 
    xlab("Significantly differentially expressed?") +
    scale_y_continuous(limits= c(0,0.0000035)) + 
    stat_pvalue_manual(data= wilcox_pwmScore,
                       label= "p.adj")
  plot(p)
  
  return(DT_list)
}


