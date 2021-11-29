# FILE:    merge_SE.R
# DATE:    5 Match 2021
# AUTHOR:  Karsten Suhre
# PURPOSE: merge to SEs
# MODIF:   9 March 2021 - remove merge_SE and generate rmerge_SE and cmerge_SE
#                         (like rbind and cbind)

rmerge_SE <- function(D1, D2) {
  
  cat("merging two SE by row:\n")
  cat("D1:\n")
  print(D1)
  cat("D2:\n")
  print(D2)
  
  # check if the colnames match
  # a future version could use left_join and pass through its options
  if (! sum(colnames(D1) == colnames(D2)) == dim(D1)[2]) {
    stop("colnames(D1) != colnames(D2)")
  }
  
  # remove all colData that does not have the same variable name and class and same content for each sample
  # if data contains NA's, it will be removed ... this should be improved in a future version
  # but I knew of no quick way of comparing vectors for identity if they contain NA's
  # (there should be a simple routine for this somewhere),
  # sum(cdata1 == cdata2, na.rm = TRUE) != dim(D1)[2] is not very elegant
  # a future version could have the option to keep non-matching variables and rename them 
  # to .1 and .2 and fill the other data with NA's
  # if the two SE's were read by Maplet, then also the merging of the meta information needs to be treated
  # i.e. for now it reports only the first loaded file in the HTML output
  for (i in intersect(names(colData(D1)),names(colData(D2))) ) {
    # cat("working on", i, "\n")
    cdata1 = D1[[i]]
    cdata2 = D2[[i]]
    lidentical = TRUE
    if (class(cdata1) != class(cdata2)) {
      lidentical = FALSE
      # cat("not identical - different class type\n")
    } else if (sum(cdata1 == cdata2, na.rm = TRUE) != dim(D1)[2]) {
      lidentical = FALSE
      # cat("not identical - different content\n")
    }
    if (! lidentical) {
      cat("deleting variable with non-matching values in D1 and D2:", i, "\n")
      colData(D1)[which(names(colData(D1)) == i)] = NULL
      colData(D2)[which(names(colData(D2)) == i)] = NULL
    }
    
  }
  
  # remove variables that are only in one data set, 
  # a future version could keep them and fill in NA's
  for (i in setdiff(names(colData(D1)),names(colData(D2))) ) {
    cat("deleting variable that is only in D1:", i, "\n")
    colData(D1)[which(names(colData(D1)) == i)] = NULL
  }  
  for (i in setdiff(names(colData(D2)),names(colData(D1))) ) {
    cat("deleting variable that is only in D2:", i, "\n")
    colData(D2)[which(names(colData(D2)) == i)] = NULL
  }  
  
  
  # treat row data: remove variables that are in only one data set
  # future: could be filled with NA's
  # also make sure they are of the same class
  
  for (i in setdiff(names(rowData(D1)),names(rowData(D2))) ) {
    cat("deleting row variable that is only in D1:", i, "\n")
    rowData(D1)[which(names(rowData(D1)) == i)] = NULL
  }
  
  for (i in setdiff(names(rowData(D2)),names(rowData(D1))) ) {
    cat("deleting row variable that is only in D2:", i, "\n")
    rowData(D2)[which(names(rowData(D2)) == i)] = NULL
  }  
  
  for (i in intersect(names(rowData(D1)),names(rowData(D2))) ) {
    # cat("working on", i, "\n")
    if (class(D1[[i]]) != class(D1[[i]])) {
      cat("deleting row variable with different class in D1 and D2:", i, "\n")
      rowData(D1)[which(names(rowData(D1)) == i)] = NULL
      rowData(D2)[which(names(rowData(D2)) == i)] = NULL
    }  
  }  
  
  # add a variable to id datasets
  rowData(D1)$original_dataset = rep("D1", dim(D1)[1])
  rowData(D2)$original_dataset = rep("D2", dim(D2)[1])
  
  # return merged SE
  rbind(D1,D2)
  
}
 

# copied the above code and edited 
cmerge_SE <- function(D1, D2) {
  
  # test code
  # D1 = D_muscle
  # D2 = D_fat
  # cmerge_SE(D_muscle, D_fat)
  
  cat("merging two SE by column:\n")
  cat("D1:\n")
  print(D1)
  cat("D2:\n")
  print(D2)
  
  # check if the colnames match
  # a future version could use left_join and pass through its options
  if (! sum(rownames(D1) == rownames(D2)) == dim(D1)[1]) {
    stop("rownames(D1) != rownames(D2)")
  }
  
  # remove all rowData that does not have the same variable name and class and same content for each sample
  for (i in intersect(names(rowData(D1)),names(rowData(D2))) ) {
    # cat("working on", i, "\n")
    cdata1 = rowData(D1)[[i]]
    cdata2 = rowData(D2)[[i]]
    lidentical = TRUE
    if (class(cdata1) != class(cdata2)) {
      lidentical = FALSE
      # cat("not identical - different class type\n")
    } else if (sum(cdata1 == cdata2, na.rm = TRUE) != dim(D1)[1]) {
      lidentical = FALSE
      # cat("not identical - different content\n")
    }
    if (! lidentical) {
      cat("deleting variable with non-matching values in D1 and D2:", i, "\n")
      rowData(D1)[which(names(rowData(D1)) == i)] = NULL
      rowData(D2)[which(names(rowData(D2)) == i)] = NULL
    }
    
  }
  
  # remove variables that are only in one data set, 
  # a future version could keep them and fill in NA's
  for (i in setdiff(names(rowData(D1)),names(rowData(D2))) ) {
    cat("deleting variable that is only in D1:", i, "\n")
    rowData(D1)[which(names(rowData(D1)) == i)] = NULL
  }  
  for (i in setdiff(names(rowData(D2)),names(rowData(D1))) ) {
    cat("deleting variable that is only in D2:", i, "\n")
    rowData(D2)[which(names(rowData(D2)) == i)] = NULL
  }  
  
  # treat row data: remove variables that are in only one data set
  # future: could be filled with NA's
  # also make sure they are of the same class
  
  for (i in setdiff(names(colData(D1)),names(colData(D2))) ) {
    cat("deleting col variable that is only in D1:", i, "\n")
    colData(D1)[which(names(colData(D1)) == i)] = NULL
  }
  
  for (i in setdiff(names(colData(D2)),names(colData(D1))) ) {
    cat("deleting col variable that is only in D2:", i, "\n")
    colData(D2)[which(names(colData(D2)) == i)] = NULL
  }  
  
  for (i in intersect(names(colData(D1)),names(colData(D2))) ) {
    # cat("working on", i, "\n")
    if (class(colData(D1)[[i]]) != class(colData(D2)[[i]])) {
      cat("deleting col variable with different class in D1 and D2:", i, "\n")
      colData(D1)[which(names(colData(D1)) == i)] = NULL
      colData(D2)[which(names(colData(D2)) == i)] = NULL
    }  
  }  
  
  # add a variable to id datasets
  colData(D1)$original_dataset = rep("D1", dim(D1)[2])
  colData(D2)$original_dataset = rep("D2", dim(D2)[2])
  
  # return merged SE
  cbind(D1,D2)
  
}
