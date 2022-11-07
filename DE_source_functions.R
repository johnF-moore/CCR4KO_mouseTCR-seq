library(edgeR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

#### Functions for DE analysis and filtering manipulation ####



DE_wrapper <- function(data_path, cell_type, counts= NULL, log2_count= NULL, overlap= NULL, logFC_threshold,
                       aa_len_min= 7, aa_len_max= 17, filt_method= NULL, save= F, out_name= NULL, savage){
  ##### Filtering ####
  dat <- read.csv(data_path) %>% 
    dplyr::select(X, contains(cell_type))
  ## Where X is the CDR3 alpha aa sequencess 
  ## and each column of reads contains the cell type in the sample name/column name
  dat <- dat[nchar(dat$X) > aa_len_min & nchar(dat$X) < aa_len_max,]
  ## Filtering for sequences that make biological sense.
  ## This filter is borrowed from Savage Lab. I don't see a reason to fiddle with it.
  cpm_matrix <- as.data.frame(cpm(y= dat[-1]))
  X <- dat$X ## Used later, but has to be pulled out before filtering.
  counts <- counts/exp(mean(log2(colSums(cpm_matrix))))
  ## normalized filtering threshold,
  ## not super necessary as long as library sizes are approximately the same for DE groups
  ## I could do RPKM normalization to correct for differences in reads due to gene length, but 
  ## for DE purposes, this isn't super necessary. Although, it has been shown to affect DE a little. 
  if(filt_method == "phenotypewise"){
    dat <- dat[(rowSums(dplyr::select(cpm_matrix, contains("CCR4KO")) >= counts) >= overlap |
                rowSums(dplyr::select(cpm_matrix, contains("WT"))     >= counts) >= overlap),]
    ## CDR3 aa will always be 1st column b/c of dplyr::select above.
    ## Hard coding "WT" and "CCR4KO"
  } else if(filt_method == "IMPACC") {
    dat <- dat[rowSums(dat[-1] > 0) >= overlap,]
  } else if(filt_method == "AvgLogCPM"){
    dat <- dat[aveLogCPM(y= dat[-1]) > log2_count &
                 rowSums(dat[-1] > 0) >= overlap,]
    ## Checking if each TCR has an average Log CPM above some threshold
    ## Could be biased by a few samples that highly express the TCR
    ## The glmQLFit should correct for that.
    ## This is the starting method employed by Aaron Lun from edgeR.
    ## Benefit is that the filtering method is independent of the experimental set up,
    ## which is a necessity and makes this filtering generalizable to multiple pipelines.
    
    ## I am adding a slight parameter that the TCR has to be expressed by more than 1 sample
    ## b/c weirdly some of these will be called DE. 
  } else if(filt_method == "by_cellType"){
    dat <- dat[rowSums(cpm_matrix >= counts) >= overlap,]
    ## Filtering by expression in the cell type rather than by condition.
  }
  
  overlap_df <- data.frame(X=dat$X,
                           num_cellType= rowSums(dat[-1] > 0)) 
  #### DE Data Prep ####
  group <- substr(x= colnames(dat[-1]), start= 1, stop= nchar(colnames(dat[-1])) -1) %>%
    factor()
  y <- DGEList(dat[,-1], group=group,
               genes=dat[,1,drop=FALSE]) 
  options(digits=3) 
  
  ## Making design matrix
  design <- model.matrix(~0+group)
  colnames(design) <- gsub(cell_type, "", levels(group))
  AveLogCPM <- aveLogCPM(y)
  hist(AveLogCPM)
  
  y <- calcNormFactors(y)
  
  if(savage){
    scaleFactor=2^mean(log2(colSums(y$counts)))
    ds=(y$counts %*% diag(scaleFactor/(y$samples$norm.factors*y$samples$lib.size)))
    colnames(ds)=colnames(y$counts)
    
    
    pseudo=min(ds[ds>0])
    
    y=estimateGLMCommonDisp(y, design)
    y=estimateGLMTrendedDisp(y,design)
    y=estimateGLMTagwiseDisp(y,design)
    
    ## Plot Data 
    pch <- c(0,1,2,15,16,17)
    colors <- rep(c("darkgreen", "red", "blue"), 2)
    plotMDS(y, col=colors[group], pch= pch[group])
    legend("topleft", legend= levels(group), pch=pch, col= colors, ncol= 2)
    
    ## Differential Expression  
    fit <- glmFit(y, design= design)
    con <- makeContrasts(CCR4KO_v_WT= CCR4KO-WT,
                               levels= design)
    tr<- glmLRT(fit, contrast= con) 
  } else{
    pch <- c(0,1,2,15,16,17)
    colors <- rep(c("darkgreen", "red", "blue"), 2)
    plotMDS(y, col=colors[group], pch=pch[group])
    legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)
    
    y <- estimateDisp(y, design, robust=TRUE, trend.method= "locfit.mixed")
    fit <- glmQLFit(y, design, robust=F, abundance.trend = F)
    #head(fit$coefficients)
    
    #### DE Test ####
    con <- makeContrasts(CCR4KO_v_WT= CCR4KO-WT, 
                         levels= design)
    tr <- glmTreat(fit, contrast= con, lfc= logFC_threshold, null = "interval") 
      ## Test logFC & FDR simultaneously
  }
  
  ## Combining data back together for output data frame
  df_tr <- as.data.frame(topTags(tr, n= nrow(tr))) 
  df_tr <- left_join(df_tr, overlap_df, by= "X")
  colnames(cpm_matrix) <- paste0(colnames(cpm_matrix),"_freq")
  cpm_matrix <- cbind(X, as.data.frame(as.matrix(cpm_matrix)/1e6))
  df_tr <- left_join(df_tr, cpm_matrix, by= "X")
  colnames(dat)[-1] <- paste0(colnames(dat)[-1], "_raw")
  df_tr <- left_join(df_tr, dat, by= "X")
  
  p <- ggplot(df_tr, aes(x= logFC, y= FDR, color= num_cellType)) + 
    geom_point(size= 2) +
    geom_hline(yintercept= 0.05,
               linetype= "dashed") + 
    geom_vline(xintercept= c(logFC_threshold, -logFC_threshold),
               linetype= "dashed") +
    scale_y_continuous(
      breaks= c(10^-6, 10^-4, 10^-2, 1),
      labels= paste(c(10^-6, 10^-4, 10^-2, 1)),
      trans= "log10",
      limits= c(10^-6, 1)) +
    scale_x_continuous() +
    xlab("logFC \n Upregulated in WT on left, CCR4KO on right") +
    scale_color_viridis_c() +
    theme(text = element_text(size= 15))
  
  plot(p)
  
  df_tr <- df_tr[df_tr$FDR < 0.05,]
  if(save){
    ggsave(plot= p, 
           filename= paste0(out_name, "_volcanoPlot"), 
           width= 9,
           height= 9)
    write.csv(df_tr, paste0(out_name, "_df.csv"),
              row.names= F)
  } ## example of a good out_name would be "/stor/work/Ehrlich/Pablo_TCRseq/DE_data/Treg_IMPACC_log2_R10_O2"
  
  is.de <- decideTestsDGE(tr, lfc= logFC_threshold)
  print(summary(is.de))
  
  return(df_tr)
}





