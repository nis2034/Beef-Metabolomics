# FILE:    beef_metabolomics_maplet_v1.R
# DATE:    23 March 2021
# AUTHOR:  Nisha Stephan
# PURPOSE: reading Biocrates p500 data and processing it in autonomics and/or maplet
#          for maplet see here: https://github.com/krumsieklab/maplet/blob/main/README.md
# MODIF:  Setup a git project for the development on 29/11/2021
#

# 

# cleanup
getwd()
#setwd("/home/rstudio/home/autonomics/nis2034/beef\ metabolomics")
setwd("/home/rstudio/host/beef_data")
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

####################################
# # load the data
####################################

# read the two data sets from Biocrates (LC-MS and FIA-MS)

# LIVER
D_FIA <- mt_load_biocrates(file = "5023_Liver-FIA_nonorm_microM.xlsx", first_data = c(17,21))
D_LC <- mt_load_biocrates(file = "5023_Liver-LC_nonorm_microM.xlsx", first_data = c(33,21))
#-------------------------------------------------------------
cat("special treatment to use taurine from run 2\n")  
zwi = D_LC[,61:120]
D_LC = D_LC[,1:60]
ix = which(rowData(D_LC)$feature_id == "Taurine")
assay(D_LC)[ix,] = assay(zwi)[ix,]
#-------------------------------------------------------------
D_liver = rmerge_SE(D1 = D_FIA, D2 = D_LC)
rowData(D_liver)$Method = ifelse(rowData(D_liver)$original_dataset == "D1", "FIA-MS", "LC-MS")

# MUSCLE
D_FIA <- mt_load_biocrates(file = "5023_Muscle-FIA_nonorm_microM.xlsx", first_data = c(17,21))
D_LC <- mt_load_biocrates(file = "5023_Muscle-LC_nonorm_microM.xlsx", first_data = c(33,21))
#-------------------------------------------------------------
cat("special treatment to use taurine from run 2\n")  
zwi = D_LC[,61:120]
D_LC = D_LC[,1:60]
ix = which(rowData(D_LC)$feature_id == "Taurine")
assay(D_LC)[ix,] = assay(zwi)[ix,]
#-------------------------------------------------------------
D_muscle = rmerge_SE(D1 = D_FIA, D2 = D_LC)
rowData(D_muscle)$Method = ifelse(rowData(D_muscle)$original_dataset == "D1", "FIA-MS", "LC-MS")

# FAT
D_FIA <- mt_load_biocrates(file = "5023_Fat-FIA_nonorm_microM.xlsx", first_data = c(17,21))
D_LC <- mt_load_biocrates(file = "5023_Fat-LC_nonorm_microM.xlsx", first_data = c(17,21))
D_fat = rmerge_SE(D1 = D_FIA, D2 = D_LC)
rowData(D_fat)$Method = ifelse(rowData(D_fat)$original_dataset == "D1", "FIA-MS", "LC-MS")


# rename the colnames and sort the D's
colnames(D_liver) = D_liver$`Tube label`

ix = sort(colnames(D_liver), index.return = TRUE)$ix 
D_liver = D_liver[,ix]

colnames(D_muscle) = D_muscle$`Tube label`

ix = sort(colnames(D_muscle), index.return = TRUE)$ix 
D_muscle = D_muscle[,ix]

colnames(D_fat) = D_fat$`Tube label`

ix = sort(colnames(D_fat), index.return = TRUE)$ix 
D_fat = D_fat[,ix]




D = cmerge_SE(D_muscle, D_fat)
zwi = ifelse(colData(D)$original_dataset == "D1", "muscle", "fat") 
D = cmerge_SE(D_liver, D)
colData(D)$tissue = ifelse(colData(D)$original_dataset == "D1", "liver", "zwi") 
colData(D)$tissue[which(colData(D)$tissue == "zwi")] = zwi

rowData(D)$tissue = rep("multiple", dim(D)[1])



# reformat the label for easier parsing (and correct some inconsistency)
helper = sub("^CH","C", colnames(D))
helper = sub("^(.)([123])","\\1.\\2", helper)
helper = sub("^(.)\\.","\\1 buf", helper)
helper = sub("buf2.1", "buf1.2", helper)
helper = sub("buf3.1", "buf1.3", helper)
helper = sub("buf1.", "Y", helper)
print(helper)

D$helper = helper

D$Material = rep("matrix", length(helper))
D$Material[grep("Y", helper)] = "buffer"
print(table(D$Material))

# Buffer described by the first character
# R = 80%Methanol + 20% PBS 
# M = 80% Methanol + 20% H2O 
# E = 50% Methanol + 50% PBS 
# CH = 50% Methanol + 50% H2O 
D$Buffer = substr(helper,1,1)
D$Buffer = sub("^R", "80% Meth + 20% PBS", D$Buffer)
D$Buffer = sub("^M", "80% Meth + 20% H2O", D$Buffer)
D$Buffer = sub("^E", "50% Meth + 50% PBS", D$Buffer)
D$Buffer = sub("^C", "50% Meth + 50% H2O", D$Buffer)
print(table(D$Buffer))

# Buffer volume description by second character:  
# A = 120 uL 
# B = 240 uL 
D$BufferVolume = substr(helper,3,3)
D$BufferVolume = sub("B", "240 uL", D$BufferVolume)
D$BufferVolume = sub("A", "120 uL", D$BufferVolume)
D$BufferVolume = sub("Y", "buffer", D$BufferVolume)
print(table(D$BufferVolume))

# get the first letter of the tube label as grouping variable
#D$Group = substr(D$`Tube label`,1,1)
D$Group = substr(helper,1,1)
print(table(D$Group))

# Replicate number described by first number after the character describing volume  
# e.g. A1 = 120uL replicate 1; B3 = 240 uL replicate 3  
D$Replicate = substr(helper,4,4)

# Tissue described by the last character 
# F = Fat; we also have L = Liver and B = muscle 
# Experiment number described by the number after character defining tissue 
# F.1 = tissue fat experiment 1; F.2 = tissue fat experiment 2 
D$Experiment = substr(helper,6,8)
D$Experiment[which(D$Experiment == "")]    = "no tissue"
D$Experiment[which(D$Experiment == "X.1")] = "exp 1"
D$Experiment[which(D$Experiment == "X.2")] = "exp 2"
print(table(D$Experiment))

print(table(D$Buffer,D$BufferVolume))
print(table(D$Buffer,D$Experiment))
print(table(D$Experiment,D$BufferVolume))
print(table(D$Experiment,D$BufferVolume, D$Buffer))
print(table(D$tissue, D$Buffer))
print(table(D$tissue, D$Experiment))

print(table(D$tissue, D$Material))

# make sure colnames are unique 
colnames(D) = paste(colnames(D)," [", colData(D)$tissue, "]", sep="")

# There cannot be 0 or negative values while log2 transforming the data 
# hence converting 0 to NA and imputing with min values
assay(D) [ which(assay(D) == 0, arr.ind = TRUE) ] = NA

D %>%  mt_reporting_data() 

# PART 1 - PREPROCESSING AND PLOTTING ----------------------------------------------------

D3 <-D %>% mt_reporting_heading(heading = "Data-Cleanup", lvl = 1) %>%  
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

# PART 2 - GET PATHWAY ANNOTATIONS ----------------------------------------------------

# rowData(D)
# 
# Dp <- D %>%
#   # heading for html file
#   mt_reporting_heading(heading = "Get Pathway Annotations", lvl = 1) %>%
#  
#   # get pathway annotations
#   #   alternative functions: mt_anno_pathways_xls, mt_anno_pathways_graphite, mt_anno_pathways_uniprot
#   mt_biocrates_anno_pathways_xls(file = "meta_anno_hmdb_h.txt", in_col = "Class", out_col = "pathway",
#                                  in_feat_id_col = "feature_id",hmdb_feat_id_col = "HMDB_id", pw_id_col = "SMP",
#                                  pw_name_col = "pathway_name",hmdb_class ="class",
#                                  hmdb_direct_parent = "direct_parent",hmdb_super_class = "super_class",
#                                  raw_db_outfile = "pathway_annotation.xlsx") %>%
#   # remove redundant
#   mt_anno_pathways_remove_redundant(feat_col = "KEGG_ids", pw_col = "pathway") %>%
#   # write pathway annotations
#   mt_write_pathways(file="ExamplePipeline_PathwayAnnotations.xlsx", pw_col = "pathway") %>%
#   {.}
# openxlsx::write.xlsx(rowData(D),"anno_biocrates.xlsx" )

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

# PART 4 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: MISSINGNESS ANALYSIS ---------------------------------------




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


