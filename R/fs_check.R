#' @export
cq_load_fcs <- function(files, is_h5 = TRUE, ...){
  res <- sapply(files, function(file)load_cytoframe_from_fcs(file, is_h5 = is_h5, ...))
  attr(res, "class") <- "cq_data"
  res
}

#' @export
print.cq_data <- function(x){
  cat("cytoqc data: \n")
  cat(length(x), " samples \n")
}

#' @export
cq_find_reference_channels <- function(cq_data, delimiter = "|"){
  keys <- sapply(cq_data, function(cf){
    paste(sort(colnames(cf)), collapse = delimiter)
  })
  res <- table(keys)
  res <- names(which.max(res))
  strsplit(res, split= delimiter, fixed = "TRUE")[[1]]
}
#' QA processes of cellcount
#'
#' detect aberations in the number of events per sample.
#' If there is a natural grouping among the samples, this can be
#' specified using the \code{grouping} argument. In this case, the
#' outlier detection will be performed within its respective group for a
#' particular sample.
#'
#' @examples
#' data(GvHD)
#' fs <- GvHD[1:10]
#' res <- fs_qc_cellcount(fs)
#' @import flowCore
#' @return a s3 object that contains the sample names, cell count and outliers detected.
#' This object can be printed or plotted by its dedicated s3 method.
fs_qc_cellcount <- function(fs, grouping=NULL, func = qoutlier, ...){

  if(!is.null(grouping))
    if(!is.character(grouping) || ! grouping %in% colnames(pData(fs)))
      stop("'grouping' must be a character indicating one of the ",
           "phenotypic variables in 'fs'")

  keyword(fs, "$TOT")

}
