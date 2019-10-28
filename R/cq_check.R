#' @export
cqc_load_fcs <- function(files, is_h5 = TRUE, ...){
  res <- sapply(files, function(file)load_cytoframe_from_fcs(file, is_h5 = is_h5, ...))
  names(res) <- basename(names(res))
  attr(res, "class") <- "cqc_data"
  res
}

#' @export
print.cqc_data <- function(x){
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
cqc_params_find_reference <- function(cqc_data, delimiter = "|", ...){
  keys <- sapply(cqc_data, function(cf){
    params <- cf_get_params(cf, ...)
    paste(sort(params), collapse = delimiter)
  })
  res <- table(keys)
  res <- names(which.max(res))
  strsplit(res, split= delimiter, fixed = "TRUE")[[1]]
}

#' @importFrom dplyr bind_rows
#' @export
cqc_params_check <- function(cqc_data, reference_params, ...){
  res <- sapply(cqc_data, function(cf){
                params <- cf_get_params(cf, ...)

                unknown <- setdiff(params, reference_params)
                missing <- setdiff(reference_params, params)
                if(length(unknown) > 0 || length(missing) > 0)
                {
                  list(unknown = unknown, missing = missing)
                }
                    }, simplify = FALSE)
  res <- Filter(Negate(is.null), res)

  attr(res, "class") <- "cqc_param_report"
  res
}
#' @export
as.data.frame.cqc_param_report <- function(x){
  res <- sapply(x, function(i){

      tibble("Not in reference" = paste(i[["unknown"]], collapse = ",")
                 , "Missing" = paste(i[["missing"]], collapse = ",")
                 )

  }, simplify = FALSE)
  bind_rows(res, .id = "FCS")
}
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @export
format.cqc_param_report <- function(x){
  if(length(x) == 0)
    x <- data.frame(FCS = "All passed")

  x <- kable(as.data.frame(x)) %>%
    kable_styling("bordered", full_width = F, position = "left") %>%
      column_spec(1, bold = TRUE) %>%
        row_spec(0, background = "#9ebcda", color = "black")
  if(length(x) == 0)
    x <- x %>% row_spec(1, color = "green")
  x
}
#' @param max.distance Maximum distance allowed for a match. See ?agrep
#' @importFrom tibble tibble add_row
#' @export
cqc_params_propose_solution <- function(check_results, max.distance = 0.1){
  res <- tibble(FCS = character(), from = character(), to = character())
  for(sn in names(check_results))
  {
    check_result <- check_results[[sn]]
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]

    #iteratively try to find the aproximate match between two vecs
    while(length(unknown) >0 && length(missing) >0)
    {
      #Levenshtein (edit) distance
      dist_mat <- adist(unknown, missing, ignore.case = TRUE)
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
      x <- unknown[ridx]
      y <- missing[cidx]
      #check if exceeds max.distance
      #agrep can be avoided if the formula of max.distance used by agrep is figured out
      ind <- agrep(x, y, ignore.case = TRUE, max.distance = max.distance)
      if(length(ind) > 0)#
      {
        #pop the matched item
        unknown <- unknown[-ridx]
        missing <- missing[-cidx]
        res <- add_row(res, FCS = sn, from = x, to = y)
      }else
        break #otherwise stop the maatching process

    }

  }

  attr(res, "class") <- c("cqc_param_solution", attr(res, "class"))
  res
}

#' @export
print.cqc_param_solution <- function(x){
  attr(x, "class") <- attr(x, "class")[-1]
  print(x)
}
#' @importFrom kableExtra collapse_rows
#' @importFrom dplyr %>% select distinct
#' @importFrom tidyr unite
#' @export
format.cqc_param_solution <- function(x, itemize = FALSE){
  if(!itemize)
    x <- x %>% select(-1) %>% distinct()

  x <- x %>%
        unite("Proposed change", from, to, sep = " --> ") %>%
          kable() %>%
            kable_styling("bordered", full_width = F, position = "left") %>%
              collapse_rows(columns = 1, "top")%>%
                row_spec(0, background = "#e5f5e0", color = "black")
  if(itemize)
    x <- x %>% column_spec(1, bold = TRUE)
  x
}

# summary.cqc_param_solution <- function(x){
#   distinct(x[,-1])
# }
#' @export
cqc_params_fix <- function(cqc_data, solution){
  invisible(apply(solution, 1, function(row){
                    sn <- row[["FCS"]]
                    cf <- cqc_data[[sn]]
                    flowWorkspace:::setChannel(cf@pointer, row[["from"]], row[["to"]])
                })
              )
}

#' @export
cqc_params_drop_not_in_reference <- function(cqc_data, check_results){
  for(sn in names(check_results))
  {
    check_result <- check_results[[sn]]
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]
    if(length(unknown) > 0 && length(missing) == 0)
    {
      cf <- cqc_data[[sn]]
      cqc_data[[sn]] <- cf[, !colnames(cf) %in% unknown]
    }
  }
  cqc_data
}

#' @export
cqc_drop_samples <- function(cqc_data, check_results){

    cqc_data[-match(names(check_results), names(cqc_data))]
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
