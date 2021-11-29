#' Add pathway information
#'
#' Adds custom pathways to the already existing \code{SummarizedExperiment} data structure using a flat file.
#' NOTE: The flat file should be an Excel file, containing 3 columns:
#' \itemize{

#'   \item{pw_id_col} contains pathway IDs
#'   \item{pw_name_col} contains pathway names}
#' The pathway columns pw_id_col and pw_name_col should be in long format, meaning if a feature is assigned to multiple pathways,
#' each feature pathway pair should appear in a separate row.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param file Path where the pathway annotation flat file is stored.
#' @param in_col rowData column to use for pathway fetching. The selected column must contain feature
#'    identifiers (e.g. HMBD, KEGG, ChEBI, etc).
#' @param out_col New column name for rowData to output pathway information to.
#' @param sheet Sheet name or number to read in flat file.

#' @param pw_id_col Name of flat file colname containing pathway IDs.
#' @param pw_name_col Name of flat file colname containing pathway names.
#' @param raw_db_outfile OPTIONAL. Name of file to export the pathway database to. Must be a string containing the path name with a
#'    .xlsx extension.
#'
#' @return rowData: New pathway annotation column added.
#' @return $results$pathways: A dataframe of pathway information.
#'
#' @examples
#' \dontrun{# annotate features using the SMP column of the pathway database flat file
#' ... %>%
#'       mt_anno_pathways_xls(file = system.file("extdata", "precalc/hmdb/hmdb_preprocessed_4.0.xlsx", package = "maplet"),
#'                            in_col = "HMDb_ID",
#'                            out_col = "janpw",
#'                            pw_id_col = "SMP",
#'                            pw_name_col = "pathway_name") %>%
#' ...}
#'
#' @author NS
#'
#' @importFrom data.table :=
#'
#' @export
require(XML)
mt_biocrates_anno_pathways_xls <- function(D,
                                 file,
                                 in_col,
                                 out_col,
                                 sheet,
                                 feat_id_col,
                                 pw_id_col,
                                 pw_name_col,
                                 raw_db_outfile) {
  
  # check arguments
  stopifnot("SummarizedExperiment" %in% class(D))
  
  if (missing(file))
    stop("file must be given to fetch pathway annotation file from")
  
  if (!file.exists(file))
    stop(glue::glue("{file} does not exist. input a valid flat file path."))
  
  if(!in_col %in% names(rowData(D)))
    stop(glue::glue("in_col is invalid. please enter an existing column name."))
  
  
  maplet:::mti_logstatus(glue::glue("reading annotation file: {basename(file)}, sheet: {sheet}"))
  # check file colnames
  pwdb <- data.frame(fread(file))
  head(pwdb)
  str(pwdb)
  pwdb$accession
  pwdb[pwdb$accession=="HMDB0000432",]
  rowData(D)
  
  hmdb_feat_id_col = "HMDB_id"
  flatfile_cols <- c(hmdb_feat_id_col, pw_id_col, pw_name_col,hmdb_class, hmdb_direct_parent,hmdb_super_class )
  valid_names <- flatfile_cols %in% names(pwdb)
  if (!all(valid_names)) {
    invalid_names <- flatfile_cols[!valid_names]
    stop(sprintf("Non-existent flat file column names. Please replace: %s \n with one of: %s",
                 stringr::str_c(invalid_names, collapse = ", "),
                 stringr::str_c(names(pwdb), collapse = ", ")))
  }
  
  pwdb %<>%
    dplyr::select(hmdb_feat_id_col = !!hmdb_feat_id_col,
                  ID = !!pw_id_col,
                  pathway_name = !!pw_name_col,
                  hmdb_class = !!hmdb_class,
                  hmdb_direct_parent = !!hmdb_direct_parent,
                  hmdb_super_class = !!hmdb_super_class)
                  
  
  # create a dataframe that enables the mapping between pathway
  # names and IDs. Included also are num_total, num_measured,
  # num_pw_total, and num_pw_measured (see below for further details)
  # in this dataframe map
  
  # num_total - overall number of features in entire database (this will
  # be a redundant, repeating number, identical in every row… but it’s the
  # easiest way to store it right now)
  # pwdb$accession is used here, since there is a many-to-one mapping
  # between HMDB_id (secondary accessions) and accession
  
  num_total <-
    pwdb$hmdb_feat_id_col %>%
    unique() %>%
    length()
  
  # num_measured - overall number of measured features found in the DB
  num_measured =
    
    # for this intersect, similar to the one in pwdb_summary (below),
    # it is assumed that the IDs in the dataset D is non-redundant.
    # else this number may not be accurate.
    rowData(D)[[in_feat_id_col]] %>% unique() %>%
    length()
  
  maplet:::mti_logstatus(glue::glue("summarizing pathway information"))
  pwdb_summary <- pwdb
  # using methods from data.table reduces runtime by almost 10x as compared
  # to dplyr
  pwdb_summary <-
    data.table::setDT(pwdb_summary)[, `:=`(
      
      num_total = num_total,
      num_measured = num_measured,
      
      # num_pw_total - the total number of features in that pathway
      # here, accession is used.
      # (overall DB background)
      num_pw_total = data.table::uniqueN(hmdb_feat_id_col),
      
      # num_pw_measured - the number of measured features in that pathway
      # (for the M type of analysis on the actual measured background).
      num_pw_measured =
        sum(hmdb_direct_parent %in% rowData(D)[[in_col]], na.rm = TRUE)
    ),
    by = ID] %>%
    # some IDs might have more than two names, however, these will be discarded
    # for now
    unique(by = c("ID")) %>%
    subset(!is.na(ID),
           select = -hmdb_feat_id_col)
  
  pwdb_summary$ID %>% unique()%>% length()
  
  maplet:::mti_logstatus(glue::glue("nesting pathway IDs"))
  # nest all the pathway IDs given our lieblings input IDs
  pwdb_reduced <-
    pwdb %>% distinct(ID, pathway_name) %>% 
    dplyr::filter(!is.na(ID)) 
  
  # match the nested pathways to our lieblings IDs
  match_idx <-
    match(rowData(D)[[in_col]],
          pwdb_reduced$pathway_name)
  
  pw_col <- pwdb_reduced$ID[match_idx]
  
  # add the pathway IDs into D
  rowData(D)[[out_col]] <- pw_col
  
  # add pathway map to the metadata of D
  metadata(D)$pathways[[out_col]] <-
    pwdb_summary
  
  
  if(!missing(raw_db_outfile)) {
    openxlsx::write.xlsx(pwdb, raw_db_outfile)
  }
  
  
  funargs <- mti_funargs()
  metadata(D)$results %<>%
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf('added pathway annotations using %s', basename(file))
    )
  
  D
}

########################### HMDB ############################################

rowData(D)[["Class"]] %>% unique()
pwdb$hmdb_direct_parent %>% unique() %>% sort()

com = pwdb$hmdb_direct_parent %>% unique() %>% intersect(rowData(D)[["Class"]] %>% unique()) %>% length()

num_measured =
  pwdb$feat_id_col %>%
 
  intersect(rowData(D)[[in_col]]) %>%
  length()
lk1 = NULL
st = NULL
library(hmdbQuery)
lk1 = HmdbEntry(prefix = "https://www.hmdb.ca/metabolites/", 
                id = "HMDB0249508")
data(hmdb1)
s= store(hmdb1)
names(s)

st = store(lk1)
names(st)

st$name
st$kegg_id
st$taxonomy


