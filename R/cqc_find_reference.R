#----------------------------------
# Functions to automatically discover the
# reference from the data by the majority count
# (to be deprecated)
#-------------------------------------------


#' @export
cqc_find_reference_channel <- function(cqc_cf_list, ...) {
  cqc_find_params_reference(cqc_cf_list, type = "channel", ...)
}

#' @export
cqc_find_reference_marker <- function(cqc_cf_list, ...) {
  cqc_find_params_reference(cqc_cf_list, type = "marker", ...)
}
#' @importFrom dplyr filter arrange
#' @importFrom tidyr separate
cqc_find_params_reference <- function(cqc_cf_list, type = c("channel", "marker"), delimiter = "|") {
  sep <- paste0(delimiter, delimiter) # double delimiter for sep params and single delimiter for sep channel and marker
  keys <- sapply(cqc_cf_list, function(cf) {
    params <- cf_get_panel(cf) %>% arrange(channel)
    if (type == "channel") {
      params <- params[["channel"]]
    } else {
      params <- filter(params, is.na(marker) == FALSE)

      params <- params[["marker"]]
    }
    paste(params, collapse = sep)
  })
  res <- table(keys)
  res <- names(which.max(res))
  res <- strsplit(res, split = sep, fixed = "TRUE")[[1]]
  res <- tibble(reference = res)
  # if(type == "marker")
  #   res <- separate(res, channel, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
  class(res) <- c("cqc_reference", class(res))
  class(res) <- c(paste("cqc_reference", type, sep = "_"), class(res))

  res
}


#' @export
cqc_remove_not_in_reference <- function(x, ...) UseMethod("cqc_remove_not_in_reference")
#' @importFrom flowWorkspace colnames
#' @export
cqc_remove_not_in_reference.cqc_report_channel <- function(x, cqc_cf_list, ...) {
  for (sn in names(x))
  {
    check_result <- x[[sn]]
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]
    if (length(unknown) > 0 && length(missing) == 0) {
      cf <- cqc_cf_list[[sn]]
      cols <- flowWorkspace::colnames(cf)
      j <- which(!cols %in% unknown)
      flowWorkspace:::subset_cytoframe_by_cols(cf@pointer, j - 1)
    }
  }
}
#' @export
cqc_remove_not_in_reference.cqc_report_marker <- function(x, cqc_cf_list, ...) {
  for (sn in names(x))
  {
    check_result <- x[[sn]]
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]
    if (length(unknown) > 0 && length(missing) == 0) {
      cf <- cqc_cf_list[[sn]]
      cf_rename_marker(cf, unknown, "")
    }
  }
}
#' @export
cqc_drop_samples <- function(cqc_cf_list, check_results) {
  cqc_cf_list[-match(names(check_results), names(cqc_cf_list))]
}
