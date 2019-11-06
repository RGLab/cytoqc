#' @export
cqc_find_reference_channel <- function(cqc_data, ...){
  cqc_find_params_reference(cqc_data, type = "channel", ...)
}

#' @export
cqc_find_reference_marker <- function(cqc_data, ...){
  cqc_find_params_reference(cqc_data, type = "marker", ...)
}
#' @importFrom dplyr filter arrange
#' @importFrom tidyr separate
cqc_find_params_reference <- function(cqc_data, type = c("channel", "marker"), delimiter = "|"){
  sep <- paste0(delimiter, delimiter)#double delimiter for sep params and single delimiter for sep channel and marker
  keys <- sapply(cqc_data, function(cf){
    params <- cf_get_params_tbl(cf) %>% arrange(channel)
    if(type == "channel")
      params <- params[["channel"]]
    else
    {
      params <- filter(params, is.na(marker) == FALSE)

      params <- params[["marker"]]

    }
    paste(params, collapse = sep)
  })
  res <- table(keys)
  res <- names(which.max(res))
  res <- strsplit(res, split= sep, fixed = "TRUE")[[1]]
  res <- tibble(reference = res)
  # if(type == "marker")
  #   res <- separate(res, channel, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
  class(res) <- c("cqc_reference", class(res))
  class(res) <- c(paste("cqc_reference", type, sep = "_"), class(res))

  res
}

