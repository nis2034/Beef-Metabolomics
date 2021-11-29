# FILE:    mt_load_olink.R
# DATE:    5 Match 2021
# AUTHOR:  Karsten Suhre
# PURPOSE: Maplet prototype reader for Biocrates mpx500 data
#
# PARAMETERS:  sheet = name of the Excel sheet with the data
#              first_data = coordinates of the Excel cell with the first assay data point
#                           presently the first data point is in row 17 and column 21 (= U17)
#                           but this may vary if Biocrates chooses to export more columns or rows,
#              format = export format, only "MetIDQ" is implemented for now, but this can be used in case Biocrates creates diverging versions of that format
#              return_only_samples = if TRUE only data for real samples (not blanks and QC samples) is returned
#
# FUTURE DEVELOPMENTS:
#
# - automatic detection of first_data is in principle possible, but maybe not needed if this remains fixed
# - parsing of the full string of cell A1, as it contains more information that may be useful
#   e.g. it is now "MetIDQ-Oxygen-DB110-3005 Concentration [µM] NoNormalization "
# - there may be obsolete or redundant variables around, like sample_id and SAMPLE_ID
# - some of the prints may be suppressed
# - add an option to select which non-sample data should be returned

mt_load_biocrates <- function (D, file, 
                               sheet = "Data Export", 
                               first_data = c(17,21),
                               format = "MetIDQ",
                               return_only_samples = TRUE)
{
  result <- list()
  if (missing(file)) 
    stop("file must be provided")
  result$info$file <- file
  if (!missing(D)) {
    if ("SummarizedExperiment" %in% class(D) == F) 
      stop("D is not of class SummarizedExperiment")
    if (length(assays(D)) != 0) 
      stop("Passed SummarizedExperiment assay must be empty!")
    result$meta <- metadata(D)
  }
  
  # compute some values
  rskip = first_data[1]-1 # rows to skip
  cskip = first_data[2]-1 # columns to skip
  
  # read the key string from the first cell to check the format
  # the full string contains more information and may be parsed in the future
  # e.g. "MetIDQ-Oxygen-DB110-3005 Concentration [µM] NoNormalization "
  key_string = readxl::read_excel(path = file, sheet = sheet, range = "A1:A1", col_names = FALSE)[[1]]
  if ( substr(key_string,1,6) != format ) { stop("input file is not in MetIDQ format")}
  cat("File format identifier is:", key_string, "\n")
  
  # read the assay data and the row data
  assay_data = readxl::read_excel(path = file, sheet = sheet, col_names = FALSE,
                                  skip = rskip, na = c("", "NA"))
  
  # split assay and row data
  row_data = assay_data[,c(1:cskip)]
  assay_data = assay_data[,-c(1:cskip)]
  
  # read the col data (which includes the row headers)
  col_data = readxl::read_excel(path = file, sheet = sheet, col_names = FALSE,
                                skip = 1, n_max =rskip-1, na = c("", "NA"), col_types ="text")
  
  # get the col header
  col_header = col_data[[cskip]]
  
  # get the row header
  row_header = unname(as.matrix(col_data[1,c(1:cskip)])[1,])
  names(row_data) = row_header
  
  # `get the col data
  col_data = col_data[,-c(1:cskip)]
  col_data = as_tibble(t(col_data))
  names(col_data) = col_header
  
  # assume that the first row contains the feature ids
  names(col_data)[1] = "feature_id"
  
  # perform a few plausibility checks and prints
  
  stopifnot(dim(row_data)[1] == dim(assay_data)[1])
  stopifnot(dim(col_data)[1] == dim(assay_data)[2])
  
  cat("Assay data for", dim(assay_data)[1], "samples and", dim(assay_data)[2], "traits read in \n")
  
  cat("There are", dim(row_data)[2], "parameters describing the samples:\n")
  print(names(row_data))
  
  cat("There are", dim(col_data)[2], "parameters describing the traits:\n")
  print(names(col_data))
  
  # create a SE & transpose the assay_data (as Maplet assumes traits to be in the columns and samples in the rows)
  assay_data = t(assay_data)
  
  sample_id_column = "Sample Bar Code"
  colnames(assay_data) = row_data[[sample_id_column]]
  
  trait_id_column = "feature_id"
  rownames(assay_data) = col_data[[trait_id_column]]
  
  D = SummarizedExperiment(assay = list(exprs = assay_data),
                           colData = row_data,
                           rowData = col_data
  )
  
  if (return_only_samples) {
    cat("Returnig only sample data\n")
    D = D[,which(D$`Sample Type` == "Sample")]
  }
  
  # create variables expected by Maplet (a few may be redundant)
  D$sample_id = colnames(D)
  D$SAMPLE_NAME = D$sample_id

  rowData(D)$name = rownames(D)
  rowData(D)$BIOCHEMICAL = rowData(D)$name
  
  cat("Returning data for", dim(D)[1], "traits and", dim(D)[2], "samples\n")
  
  if (length(which(duplicated(colnames(D)))) >0 ) {
    cat("WARNING: duplicated sample ids\n")
    print(colnames(D)[which(duplicated(colnames(D)))])
  }
  
  if (length(which(duplicated(rownames(D)))) >0 ) {
    cat("WARNING: duplicated traits ids\n")
    print(rownames(D)[which(duplicated(rownames(D)))])
  }
  
  if (!is.null(result$meta$results)) 
    metadata(D)$results <- result$meta$results
  if (!is.null(result$meta$settings)) 
    metadata(D)$settings <- result$meta$settings

  funargs <- maplet:::mti_funargs()
  metadata(D)$results %<>% maplet:::mti_generate_result(funargs = funargs, 
                                               logtxt = sprintf("loaded Biocrates file: %s", basename(file)))
  D
}
