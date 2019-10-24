#' @export
cq_load_fcs <- function(files, is_h5 = TRUE, ...){
  res <- sapply(files, function(file)load_cytoframe_from_fcs(file, is_h5 = is_h5, ...))
  names(res) <- basename(names(res))
  attr(res, "class") <- "cq_data"
  res
}

#' @export
print.cq_data <- function(x){
  cat("cytoqc data: \n")
  cat(length(x), " samples \n")
}

cf_get_params <- function(cf, type = c("channel", "marker", "both")){
  type <- match.arg(type)
  pd <- pData(parameters(cf))
  channels <- pd[["name"]]
  markers <- pd[["desc"]]
  if(type == "channel")
    params <- channels
  else if(type == "marker")
    params <- markers[!is.na(markers)]
  else
  {
    stop("type = 'both' is not supported yet!")
    markers[is.na(markers)] <- ""
    params <- paste(channels, markers, sep = ":")

  }
  params
}
#' @export
cq_find_reference_params <- function(cq_data, delimiter = "|", ...){
  keys <- sapply(cq_data, function(cf){
    params <- cf_get_params(cf, ...)
    paste(sort(params), collapse = delimiter)
  })
  res <- table(keys)
  res <- names(which.max(res))
  strsplit(res, split= delimiter, fixed = "TRUE")[[1]]
}

#' @importFrom dplyr bind_rows
#' @export
cq_check_params <- function(cq_data, reference_params, ...){
  res <- sapply(cq_data, function(cf){
                params <- cf_get_params(cf, ...)

                unmatched <- setdiff(params, reference_params)
                missing <- setdiff(reference_params, params)
                if(length(unmatched) > 0 || length(missing) > 0)
                {
                  list(unmatched = unmatched, missing = missing)
                }
                    }, simplify = FALSE)
  res <- Filter(Negate(is.null), res)

  attr(res, "class") <- "cq_param_report"
  res
}
#' @export
print.cq_param_report <- function(x){
  res <- sapply(x, function(i){

      data.frame(unmatched = paste(i[["unmatched"]], collapse = ",")
                 , missing = paste(i[["missing"]], collapse = ",")
                 , stringsAsFactors = FALSE)

  }, simplify = FALSE)
  print(bind_rows(res, .id = "FCS"))
}

#' @param max.distance Maximum distance allowed for a match. See ?agrep
#' @export
cq_fix_param_solution <- function(check_results, max.distance = 0.1){
  res <- sapply(check_results, function(check_result){

        unmatched <- check_result[["unmatched"]]
        missing <- check_result[["missing"]]
        solution <- NULL
        #iteratively try to find the aproximate match between two vecs
        while(length(unmatched) >0 && length(missing) >0)
        {
          #Levenshtein (edit) distance
          dist_mat <- adist(unmatched, missing, ignore.case = TRUE)
          nrows <- nrow(dist_mat)
          ncols <- ncol(dist_mat)
          #pick the best match
          idx <- which.min(dist_mat)
          # browser()
          #get x, y coordinates
          ridx <- idx %% nrows
          if(ridx==0)
            ridx = nrows
          cidx <- ceiling(idx / nrows)
          #get the pair
          x <- unmatched[ridx]
          y <- missing[cidx]
          #check if exceeds max.distance
          #agrep can be avoided if the formula of max.distance used by agrep is figured out
          ind <- agrep(x, y, ignore.case = TRUE, max.distance = max.distance)
          if(length(ind) > 0)#
          {
            #pop the matched item
            unmatched <- unmatched[-ridx]
            missing <- missing[-cidx]
            solution <- rbind(solution, data.frame(from = x, to = y, stringsAsFactors = FALSE))
          }else
            break #otherwise stop the maatching process

        }
        solution
      })
   res <- bind_rows(res, .id = "FCS")
  # res <- Filter(Negate(is.null), res)
  #
  # attr(res, "class") <- c(attr(res, "class"), "cq_param_solution")
  res
}

# summary.cq_param_solution <- function(x){
#   distinct(x[,-1])
# }
cq_fix_params <- function(cq_data, solution){

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
#' @import flowCore flowWorkspace
#' @return a s3 object that contains the sample names, cell count and outliers detected.
#' This object can be printed or plotted by its dedicated s3 method.
fs_qc_cellcount <- function(fs, grouping=NULL, func = qoutlier, ...){

  if(!is.null(grouping))
    if(!is.character(grouping) || ! grouping %in% colnames(pData(fs)))
      stop("'grouping' must be a character indicating one of the ",
           "phenotypic variables in 'fs'")

  keyword(fs, "$TOT")

}
