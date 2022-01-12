# # ------------ maplet Example Pipeline ------------
# FILE:    beef_metabolomics_maplet_v1.R
# DATE:    23 March 2021
# AUTHOR:  Nisha Stephan
# PURPOSE: reading Biocrates p500 data and processing it in autonomics and/or maplet
#          for maplet see here: https://github.com/krumsieklab/maplet/blob/main/README.md
# MODIF:  Setup a git project for the development on 29/11/2021
# This example pipeline demonstrates the use of 79 of the 86 functions available from maplet.
# cleanup
getwd()
setwd("/home/rstudio/host/beef_git_tensor/Beef-Metabolomics")
rm(list=ls())

# libraries
library(tidyverse)
library(autonomics)
library(maplet)
library(magrittr)
library(grid)
library(gridExtra)
library(data.table)
library(xlsx)
# source prototype functions (should appear in Maplet one day)
source("mt_read_biocrates.R")
source("merge_SE.R")
source("mt_remove_buffer.R")
source("mt_stats_univ_lm_temp.R")
source("mti_format_se_samplewise_temp.R")
source("bt_biocrates_load_preprocess.R")
source("mt_biocrates_anno_pathways_xls")

# source("mt_stats_univ_wilcox_temp.R")
# choose options here
lremove_buffer = TRUE

# choose the tissue
ctissue = "fat"
ctissue = "liver"
ctissue = "muscle"
ctissue = "rmerge" # merge all three replicating metabolites, e.g. C2 [fat], C2 [liver], ...
ctissue = "cmerge" # merge all three replicating experiments, e.g. exp1 [fat], exp1 [liver], ...

# define basename for output
basename = "beef_metabolomics_maplet_v1"
basename = paste(basename, "_", ctissue, sep = "")
if (lremove_buffer) {
  basename = paste(basename, "_buffer_removed", sep = "")
} else {
  basename = paste(basename, "_incl_buffer", sep = "")
}

# open a log file
logfile = paste(basename, ".log", sep = "")
options(warn=-1); sink(); options(warn=0)
sink(file = logfile, split = TRUE)
cat("Start logging to", logfile, "\n")
print(Sys.time())

D <- bt_biocrates_load_preprocess()

D %>%  mt_reporting_data() 

# PART 1 - PREPROCESSING AND PLOTTING ----------------------------------------------------

D <-D %>% mt_reporting_heading(heading = "Data-Cleanup", lvl = 1) %>%  
    mt_reporting_heading(heading = "Missingness", lvl = 2) %>%
    mt_reporting_text(text = "Missigness before removing buffer") %>%
    mt_plots_missingness(feat_max=0.5) %>%
    mt_remove_buffer() %>%
    mt_reporting_text(text = "Missingness after removing buffer") %>%
    mt_plots_missingness(feat_max=0.5) %>%
    # filter metabolites with more than 50% missing values per group, otherwise PCA will throw error
    mt_pre_filter_missingness(feat_max = 0.5, group_col = "Group") %>%
    mt_reporting_text(text = "Missingness after filtering metabolites with more than 50% missing values per group") %>%
    # plot missingness distribution after filtering
    mt_plots_missingness(feat_max=0.5) %>%
    mt_reporting_heading(heading = "Sample Box Plots", lvl = 2) %>%
    mt_plots_sample_boxplot(color=Buffer, title = "Buffer Types", plot_logged = T) %>%
    mt_plots_sample_boxplot(color=tissue, title = "Tissue Types", plot_logged = T) %>%
    #' Quotient normalization
    #' Implementation according to Dieterle et al., 2006
    #'Normalize abundances using probabilistic quotient
    mt_pre_norm_quot(feat_max = 0.2, ref_samples = Buffer=="80% Meth + 20% H2O") %>%
    # show dilution plot
    mt_plots_dilution_factor(in_col="Buffer") %>%
    # plot sample boxplots after normalization
    mt_plots_sample_boxplot(color=Buffer, title = "After probabilistic quotient normalization", plot_logged = T) %>%
    # There cannot be 0 in the assay
    
    # log transform -  log2 transformation
    #   other data transformation functions: mt_pre_trans_exp, mt_pre_trans_relative, mt_pre_trans_scale
    mt_pre_trans_log() %>%
  
    mt_plots_sample_boxplot(color=Buffer, title = "After log2 transformation", plot_logged = T) %>%
  
    # impute missing values using KNN method : k Number of nearest neighbors to consider. 
    #   alternative imputation functions: mt_pre_impute_min
    # mt_pre_impute_min() %>%
    # plot sample boxplot after imputation
    ##mt_plots_sample_boxplot(color=Buffer, title = "After knn imputation", plot_logged = T) %>%
    # outlier detection (univariate)
    #   alternative functions: mt_pre_outlier_detection_mahalanobis(), mt_pre_outlier_detection_leverage()
    #   related function: mt_pre_outlier_to_na()
    # 2 columns "outlier" and "score" will be added to colData(D)
    mt_pre_outlier_detection_univariate() %>%
    # print infos about dataset
    mt_reporting_data() %>%
    # write preprocessed data to Excel file
    #   other writing functions: mt_write_se_rds (save SummarizedExerpiment object)
    mt_write_se_xls(file = "PreprocessedData.xlsx") %>%
  {.}

assay(D)[ which(is.na(assay(D)), arr.ind = TRUE) ] = 0

# PART 3 - GLOBAL STATISTICS ----------------------------------------------------

D <- D %>%
  
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  # plot PCA for different buffers
  #mt_plots_pca(scale_data = T, title = "scaled PCA - Buffer", color=Buffer, size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP for different buffers
  mt_plots_umap(scale_data = T, title = "scaled UMAP - Buffer", color=Buffer, size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap buffer
  mt_plots_heatmap(scale_data = T, annotation_col = c("Buffer"), annotation_row = c("Class"),
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  # plot PCA for different tissue types
  #mt_plots_pca(scale_data = T, title = "scaled PCA", color=tissue, size=2.5, ggadd=scale_size_identity()) %>%
  # plot UMAP for different tissue types
  mt_plots_umap(scale_data = T, title = "scaled UMAP - tissue", color=tissue, size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap tissue types
  mt_plots_heatmap(scale_data = T, annotation_col = c("tissue"), annotation_row = c("Class"),
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  
  # plot pca for buffer and tissue types together in one plot
  mt_plots_pca(scale_data = T, title = "scaled PCA", color= Buffer, shape = tissue, size=2.5, ggadd=scale_size_identity()) %>%
  
  {.}


# PART 5 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: LINEAR REGRESSION (t-test) ----------------------------
rownames(D) = make.names(rownames(D))
D$factor = as.factor(D$`Sample Description`)
# Dstat <- D[which(apply(assay(D), 1, sd)==0),] 

# rownames(Dstat)
D1 <- D

D1 <- D1 %>% 
  # heading for html file
  mt_reporting_heading(heading = "Buffer analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  #   alternative functions: mt_stats_univ_wilcox, mt_stats_univ_lm_matrixeqtl
  mt_stats_univ_lm(formula = ~ Buffer,
                   stat_name = "Buffer met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "Buffer met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Buffer met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Buffer met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "BufferAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Buffer met") %>%
  # volcano plot as overview of results
  # mt_plots_volcano(stat_name = "Buffer met",
  #                  x = fc,
  #                  feat_filter = p.adj < 0.05,
  #                  colour       = p.adj < 0.05) %>%
  # boxplot
  mt_plots_box_scatter(stat_name          ="Buffer met",
                       x                  = Buffer,
                       fill               = Buffer,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  {.}


D1 <- D1 %>% 
  # heading for html file
  mt_reporting_heading(heading = "Tissue analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  #   alternative functions: mt_stats_univ_wilcox, mt_stats_univ_lm_matrixeqtl
  mt_stats_univ_lm(formula = ~ tissue,
                   stat_name = "Tissue met") %>%
  # add fold change
  #mt_post_fold_change(stat_name = "Buffer met") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Tissue met", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Tissue met", stat_filter = p.adj < 0.05) %>%
  # write statistical results to file
  mt_write_stats(file = "TissueAnalysis.xlsx") %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Tissue met") %>%
  # volcano plot as overview of results
  # mt_plots_volcano(stat_name = "Buffer met",
  #                  x = fc,
  #                  feat_filter = p.adj < 0.05,
  #                  colour       = p.adj < 0.05) %>%
  # boxplot
 
  {.}

# PART 7.1 - STATISTICAL RESULTS PRESENTATION: STATS PATHWAY BAR PLOT & PATHVIEW ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Statisitcal Results Presentation", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Barplot", lvl = 2) %>%
  # create statsbarplots
  mt_plots_stats_pathway_bar(stat_list = c("Buffer met","Tissue met"),
                             feat_filter = p.adj < 0.05,
                             group_col = "Class",
                             
                             y_scale = "count",
                             sort_by_y = T,
                             assoc_sign_col = "statistic") %>%
  mt_plots_stats_pathway_bar(stat_list = c("Buffer met","Tissue met"),
                             feat_filter = p.adj < 0.05,
                             group_col = "Class",
                            
                             y_scale = "fraction",
                             sort_by_y = T,
                             assoc_sign_col = "statistic") %>%
 
  {.}


# PART 7.2 - STATISTICAL RESULT PRESENTATION: MULTIPLE STATISTICS HEATMAP --------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(cutoff = 0.05) %>%
  {.}

# PART 7.3 - STATISTICAL RESULT PRESENTATION: RESULT COMPARISON ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Result Comparison", lvl = 2) %>%
  # comparison plot
  mt_plots_stats_compare(stat1 = "Buffer met", filter1 = p.adj < 0.05,
                         D2 = D1, stat2 = "Tissue met", filter2 = p.adj < 0.05,
                         filter_op = "OR") %>%
  {.}


# PART 8 - PARTIAL CORRELATION NETWORK ----------------------------------------------------

D1 <- D1 %>%
  # heading for html file
  mt_reporting_heading(heading = "Partial Correlation Network", lvl = 2) %>%
  # compute partial correlation matrix
  mt_stats_cormat_genenet(stat_name = "GGM") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "GGM", method = "BH") %>%
  # plot network and color according to age analysis
  # NOTE: THIS FUNCTION CAN TAKE SEVERAL MINUTES TO RUN
  mt_plots_net(stat_name = "GGM", cor_filter = p.adj < 0.05, node_coloring = "Buffer met") %>%
  {.}



# PART 9 - PATHWAY AGGREGATION ANALYSIS ----------------------------------------------------

# NOTE ON BEST PRACTICES: This is now first aggregating the metabolite matrix into pathways creating a new matrix of pathway
#   concentration values, and then repeating the parts of the same pipeline as above. In a real scenario, you would not
#   perform both kinds of analysis within the same pipeline.

# create another SE object for second analysis branch (pathways)
D2 <- D
colnames(rowData(D2)) = make.names(colnames(rowData(D2)))

D2 <- D2 %>%
  # aggregate metabolites in the same pathways
  mt_modify_agg_pathways(pw_col = "Class", method = "aggmean") %>%
  
  # STATISTICAL ANALYSIS, OUTCOME: AGE
  # heading for html file
  mt_reporting_heading(heading = "Pathway Aggregation Analysis", lvl = 1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Buffer analysis", lvl = 2) %>%
  mt_stats_univ_lm(formula = ~ Buffer,
                   stat_name = "Buffer pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Buffer pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Buffer pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Buffer pw") %>%
  # volcano plot as overview of results
  mt_plots_box_scatter(stat_name          ="Buffer pw",
                       x                  = Buffer,
                       fill               = Buffer,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") %>%
  
   mt_reporting_heading(heading = "Tissue analysis", lvl = 2) %>%
  # linear model for binary function (equivalent to t-test)
  #   alternative functions: mt_stats_univ_wilcox, mt_stats_univ_lm_matrixeqtl
  mt_stats_univ_lm(formula = ~ tissue,
                   stat_name = "Tissue pw") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "Tissue pw", method = "BH") %>%
  # add stats logging
  mt_reporting_stats(stat_name = "Tissue pw", stat_filter = p.adj < 0.05) %>%
  # pvalue histogram
  mt_plots_pval_hist(stat_list = "Tissue pw") %>%
  # volcano plot as overview of results
  mt_plots_box_scatter(stat_name          ="Tissue pw",
                       x                  = tissue,
                       fill               = tissue,
                       plot_type          = "box",
                       feat_filter       = p.adj < 0.05,
                       feat_sort         = p.value,
                       annotation         = "{sprintf('P-value: %.2e', p.value)}\nPadj: {sprintf('%.2e', p.adj)}") 
  


D2 <- D2 %>%
  # heading for html file
  mt_reporting_heading(heading = "Multiple Statistics Heatmap", lvl = 2) %>%
  # heatmap of all statistical results
  mt_plots_multstats_heatmap(cutoff = 0.05) %>%
  # heading for html file
  mt_reporting_heading(heading = "Result Comparison", lvl = 2) %>%
  # comparison plot
  mt_plots_stats_compare(stat1 = "Buffer pw", filter1 = p.adj < 0.05,
                         D2 = D2, stat2 = "Tissue pw", filter2 = p.adj < 0.05,
                         filter_op = "OR") 


################### Generate html ###############################################
# PART 13 - CREATE ANALYSIS REPORTS ----------------------------------------------------
 D1 %>% mt_reporting_html(file = paste(basename, ".maplet.html", sep = ""),
                     title = basename) 
#   
 D2 %>% mt_reporting_html(file = "Example_Pipeline_Pathway_Analysis.html",
                          title = "Example Pipeline - Pathway Aggregation Analysis")

print(Sys.time())
sink()










# #####################%%%%%%%%%%%%%%%%%%%%%%%%%%set Buffer as subgroup variable


