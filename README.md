# CCR4-knockout mouse TCR-seq
**Bulk T cell receptor (TCR) sequencing of naive mouse T cells from the thymus of CCR4-knockout and wild-type mice**
The mice have a fixed beta-chain in the TCR that keeps half of the protein constant. This reduces the possible number of TCR sequences that can be found because only some TCR alpha-chains can pair with the fixed beta-chain. 

- DE_source_functions.R and vdj_gene_comparison_functions.R contain functions that I use for all of my differential expression and TCR analyses. 

- Finding differences in V gene and J gene usage in the TCRs of wild-type and CCR4-knockout mouse is done using comparing_VJuCDR3_normalized_gene_usage.qmd

- comparing_aa_comp_bw_WT_CCR4KO.qmd compares amino acid composition by position between WT & CCR4KO TCRs (See images below.)


Amino acid usage by position for each cell type and genotype 
![Screen Shot 2022-11-07 at 12 21 55 AM](https://user-images.githubusercontent.com/98127654/200244113-d02c6132-be8d-4635-9a33-09e48c8fe951.png)

Maximum difference in amino acid proportion between CCR4KO and wild-type mice shows that some positions are more variable, even across TCRs of different lengths.
![Screen Shot 2022-11-07 at 12 53 49 AM](https://user-images.githubusercontent.com/98127654/200244058-642470f6-16ab-4d35-9caa-0c7ed5a604e9.png)
