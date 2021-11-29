mti_format_se_samplewise_temp <- function(D){
  # coldata and assay cannot have overlapping names
  # (this should be caugh earlier in the pipeline, but here is where it causes trouble)
  inters <- intersect(colnames(colData(D)), rownames(D))
  if (length(inters)>0) {
    stop(sprintf("There are metabolites and colData variables with the same name: %s", paste0(inters, collapse = ", ")))
  }
  # cbind
  cbind(colData(D),
        t(assay(D))) %>%
    data.frame(check.names=FALSE) %>%
    tibble::rownames_to_column("merge.primary")
}
