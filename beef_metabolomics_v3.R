# FILE:    beef_metabolomics_v2.R
# DATE:    4 March 2021
# AUTHOR:  Karsten Suhre
# PURPOSE: reading Biocrates p500 data and processing it in autonomics and/or maplet
#          for maplet see here: https://github.com/krumsieklab/maplet/blob/main/README.md
# MODIF:   5 Match - add read and merge functions
#                    update phenotype definition
#          9 March 2021 - add Liver and Muscle data
#
# Buffer described by the first character
# R = 80%Methanol + 20% PBS 
# M = 80% Methanol + 20% H2O 
# E = 50% Methanol + 50% PBS 
# CH = 50% Methanol + 50% H2O 
#
# Buffer volume description by second character:  
# A = 120 uL 
# B = 240 uL 
# 
# Replicate number described by first number after the character describing volume  
# e.g. A1 = 120uL replicate 1; B3 = 240 uL replicate 3  
#
# Tissue described by the last character 
# F = Fat; we also have L = Liver and B = muscle 
# Experiment number described by the number after character defining tissue 
# F.1 = tissue fat experiment 1; F.2 = tissue fat experiment 2 
#
# BUGS:
#   non UTF-8 characters in file names kill the output of mt_reporting_html
#   FIAfile = "5023_Fat-FIA_nonorm_µM.xlsx"
#   LCfile = "5023_Fat-LC_nonorm_µM.xlsx"
#   work-around - rename the file

# cleanup
rm(list=ls())

# libraries
library(tidyverse)
#library(autonomics)
library(maplet)

# source prototype functions (should appear in Maplet one day)
source("mt_read_biocrates.R")
source("merge_SE.R")

# choose options here
lremove_buffer = FALSE

# choose the tissue
ctissue = "fat"
ctissue = "liver"
ctissue = "muscle"
ctissue = "rmerge" # merge all three replicating metabolites, e.g. C2 [fat], C2 [liver], ...
ctissue = "cmerge" # merge all three replicating experiments, e.g. exp1 [fat], exp1 [liver], ...

# define basename for output
basename = "beef_metabolomics_v3"
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
if (ctissue == "rmerge") {
  colnames(D_liver) = sub("L\\.","X.",colnames(D_liver))
}
ix = sort(colnames(D_liver), index.return = TRUE)$ix 
D_liver = D_liver[,ix]

colnames(D_muscle) = D_muscle$`Tube label`
if (ctissue == "rmerge") {
  colnames(D_muscle) = sub("B\\.","X.",colnames(D_muscle) )
}
ix = sort(colnames(D_muscle), index.return = TRUE)$ix 
D_muscle = D_muscle[,ix]

colnames(D_fat) = D_fat$`Tube label`
if (ctissue == "rmerge") {
   colnames(D_fat) = sub("F\\.","X.",colnames(D_fat))
}
ix = sort(colnames(D_fat), index.return = TRUE)$ix 
D_fat = D_fat[,ix]

# select dataset
if (ctissue == "liver") {
  D = D_liver
  rowData(D)$tissue = rep("liver", dim(D)[1])
  colData(D)$tissue = rep("liver", dim(D)[2])
} else if (ctissue == "muscle") {
  D = D_muscle
  rowData(D)$tissue = rep("muscle", dim(D)[1])
  colData(D)$tissue = rep("muscle", dim(D)[2])
} else if (ctissue == "fat") {
  D = D_fat
  rowData(D)$tissue = rep("fat", dim(D)[1])
  colData(D)$tissue = rep("fat", dim(D)[2])
} else if (ctissue == "rmerge") {
  
  # create unique rownames
  rownames(D_liver) = paste(rownames(D_liver), "[liver]")
  rownames(D_muscle) = paste(rownames(D_muscle), "[muscle]")
  rownames(D_fat) = paste(rownames(D_fat), "[fat]")
  
  D = rmerge_SE(D_muscle, D_fat)
  zwi = ifelse(rowData(D)$original_dataset == "D1", "muscle", "fat") 
  D = rmerge_SE(D_liver, D)
  rowData(D)$tissue = ifelse(rowData(D)$original_dataset == "D1", "liver", "zwi") 
  rowData(D)$tissue[which(rowData(D)$tissue == "zwi")] = zwi

  colData(D)$tissue = rep("multiple", dim(D)[2])
  
} else if (ctissue == "cmerge") {
  
  D = cmerge_SE(D_muscle, D_fat)
  zwi = ifelse(colData(D)$original_dataset == "D1", "muscle", "fat") 
  D = cmerge_SE(D_liver, D)
  colData(D)$tissue = ifelse(colData(D)$original_dataset == "D1", "liver", "zwi") 
  colData(D)$tissue[which(colData(D)$tissue == "zwi")] = zwi

  rowData(D)$tissue = rep("multiple", dim(D)[1])
  
} else {
  stop("tissue '", ctissue, "' not suported")
}

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

# remove buffer
if (lremove_buffer) {
  cat("removing buffer samples\n")
  D = D[,-which(D$Material == "buffer")]
}

# set traits with zero values to missing for missing value detection 
# to remove samples with too many zeros, will be set back to zero after that step
assay(D) [ which(assay(D) == 0, arr.ind = TRUE) ] = NA

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
  # heading for html file
  mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Zero values will be considered as missing.") %>%
  # plot missingness distribution
  mt_plots_missingness(feat_max=0.5) %>%
  # filter metabolites with more than 50% missing values per group
  mt_pre_filter_missingness(feat_max = 0.5, group_col = "Group") %>%
  # plot missingness distribution after filtering
  mt_plots_missingness(feat_max=0.5) %>%
  # add missingness percentage as annotation to samples (remaining missing)
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
  {.}

# set traits with NA values back to zero
assay(D) [ which(is.na(assay(D)), arr.ind = TRUE) ] = 0

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # section text
  mt_reporting_text(text = "Plot sample boxplots, write pre-processed data to file.") %>%
  # plot sample boxplots
  mt_plots_sample_boxplot(color=Material, title = "Original", plot_logged = T) %>%
  mt_plots_sample_boxplot(color=tissue, title = "Original", plot_logged = T) %>%
  # outlier detection (univariate)
  mt_pre_outlier_detection_univariate() %>%
  # print infos about dataset
  mt_reporting_data() %>%
  # write preprocessed data to Excel file
  #   other writing functions: mt_write_se_rds (save SummarizedExerpiment object)
  mt_write_se_xls(file = paste(basename, ".preprocessed.xlsx", sep = "")) %>%
  {.}
  
D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  # plot PCA
  mt_plots_pca(scale_data = T, title = "scaled PCA", color=Material, size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_pca(scale_data = T, title = "scaled PCA", symbol=Buffer, color=tissue, size=2.5, ggadd=scale_size_identity()) %>%
  # plot heatmap
 mt_plots_heatmap(scale_data = T, annotation_col = c("tissue"), annotation_row = c("Class"),
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  {.}

# PART 8 - PARTIAL CORRELATION NETWORK ----------------------------------------------------

D <- D %>%
  # heading for html file
  mt_reporting_heading(heading = "Partial Correlation Network", lvl = 2) %>%
  # compute partial correlation matrix
  mt_stats_cormat_genenet(stat_name = "GGM") %>%
  # add multiple testing correction
  mt_post_multtest(stat_name = "GGM", method = "BH") %>%
  # plot network and color according to age analysis
  # NOTE: THIS FUNCTION CAN TAKE SEVERAL MINUTES TO RUN
  # mt_plots_net(stat_name = "GGM", cor_filter = p.adj < 0.05) %>%
  {.}


D %>% mt_reporting_html(file = paste(basename, ".maplet.html", sep = ""),
                         title = basename)


# close log file
print(Sys.time())
sink()