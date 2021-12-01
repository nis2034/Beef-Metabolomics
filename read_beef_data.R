read_beef_data <- function() {
  
  
  
  
  
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
  
  D
  
}





D <- read_beef_data()
D %<>% mt_plots_missingness(feat_max=0.5)