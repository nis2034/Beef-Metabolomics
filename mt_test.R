library(tidyverse)
library(autonomics)
library(maplet)
library(magrittr)
library(grid)
library(gridExtra)
library(data.table)
library(shiny)
library("shinydashboard")
# source prototype functions (should appear in Maplet one day)
source("mt_read_biocrates.R")
source("merge_SE.R")
source("mt_remove_buffer.R")
source("mt_stats_univ_lm_temp.R")
source("mti_format_se_samplewise_temp.R")


source("beef_metabolomics_maplet_fun.R")

source("read_beef_data.R")

D <- read_beef_data()
D %<>% mt_plots_missingness(feat_max=0.5)

r <- metadata(D)$results

r[[1]]$fun[1]
r[[1]]$output[1]
r[[1]]$output[1]
if (r[[i]]$fun[1]=="plots")