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
#' @importFrom tibble as.tibble
#' @importFrom dplyr select rename
cf_get_params_tbl <- function(cf){
  pData(parameters(cf)) %>%
    as.tibble() %>%
      select(c("name", "desc")) %>%
        rename(channel = name, marker = desc)

}

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

#' @export
format.cqc_reference <- function(x){
    kable(x) %>%
      kable_styling("bordered", full_width = F, position = "left") %>%
          row_spec(0, background = "gray", color = "black")
}


#' @export
cqc_check <- function(x, ...)UseMethod("cqc_check")

#' @export
cqc_check.cqc_reference_channel <- function(x, ...){
  res <- cqc_check_params(x, type = "channel", ...)
  class(res) <- c("cqc_report_channel", class(res))
  res
}

#' @export
cqc_check.cqc_reference_marker <- function(x, ...){
  res <- cqc_check_params(x, type = "marker", ...)
  class(res) <- c("cqc_report_marker", class(res))
  res
}

#' @importFrom dplyr bind_rows
cqc_check_params <- function(reference_params, type, cqc_data, delimiter ="|"){
  res <- sapply(cqc_data, function(cf){
                params <- cf_get_params_tbl(cf)
                if(type == "marker")
                  params <- filter(params, is.na(marker) == FALSE)

                data <- params[[type]]
                ref <- reference_params[["reference"]]
                unknown <- setdiff(data, ref)
                missing <- setdiff(ref, data)
                if(length(unknown) > 0 || length(missing) > 0)
                {
                  list(unknown = unknown, missing = missing)
                }
                    }, simplify = FALSE)
  res <- Filter(Negate(is.null), res)

  class(res) <- c("cqc_report_params", class(res))
  res
}
#' @export
as.data.frame.cqc_report_params <- function(x){
  res <- sapply(x, function(i){

      tibble("Not in reference" = paste(i[["unknown"]], collapse = ",")
                 , "Missing channels" = paste(i[["missing"]], collapse = ",")
                 )

  }, simplify = FALSE)
  bind_rows(res, .id = "FCS")
}
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @export
format.cqc_report_params <- function(x){
  if(length(x) == 0)
    res <- kable(data.frame(FCS = "All passed"), col.names = NULL)%>% row_spec(1, color = "green")
  else
    res <- kable(as.data.frame(x)) %>%
      row_spec(0, background = "#9ebcda", color = "black")

  res <- res %>%
    kable_styling("bordered", full_width = F, position = "left") %>%
      column_spec(1, bold = TRUE)
  res
}

#' @export
cqc_find_solution <- function(x, ...)UseMethod("cqc_find_solution")
#' @export
cqc_find_solution.cqc_report_channel <- function(x, ...){
  res <- cqc_find_solution.cqc_report(x, ...)
  attr(res, "class") <- c("cqc_solution_channel", attr(res, "class"))
  res

}
#' @export
cqc_find_solution.cqc_report_marker <- function(x, ...){
  res <- cqc_find_solution.cqc_report(x, ...)
  attr(res, "class") <- c("cqc_solution_marker", attr(res, "class"))
  res

}

#' @param max.distance Maximum distance allowed for a match. See ?agrep
#' @importFrom tibble tibble add_row
#' @export
cqc_find_solution.cqc_report <- function(x, max.distance = 0.1){
  res <- tibble(FCS = character(), from = character(), to = character())
  for(sn in names(x))
  {
    check_result <- x[[sn]]
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
      from <- unknown[ridx]
      to <- missing[cidx]
      #check if exceeds max.distance
      #agrep can be avoided if the formula of max.distance used by agrep is figured out
      ind <- agrep(from, to, ignore.case = TRUE, max.distance = max.distance)
      if(length(ind) > 0)#
      {
        #pop the matched item
        unknown <- unknown[-ridx]
        missing <- missing[-cidx]
        res <- add_row(res, FCS = sn, from = from, to = to)
      }else
        break #otherwise stop the maatching process

    }

  }
  attr(res, "class") <- c("cqc_solution", attr(res, "class"))
  res
}

#' @export
print.cqc_solution <- function(x){
  attr(x, "class") <- attr(x, "class")[-(1:2)]
  print(x)
}
#' @importFrom kableExtra collapse_rows
#' @importFrom dplyr %>% select distinct
#' @importFrom tidyr unite
#' @export
format.cqc_solution <- function(x, itemize = FALSE){
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

#' @export
cqc_fix <- function(x, ...)UseMethod("cqc_fix")
#' @export
cqc_fix.cqc_solution_channel <- function(x, cqc_data){
  invisible(apply(x, 1, function(row){
                    sn <- row[["FCS"]]
                    cf <- cqc_data[[sn]]
                    flowWorkspace:::setChannel(cf@pointer, row[["from"]], row[["to"]])
                })
              )
}

#' @export
cqc_fix.cqc_solution_marker <- function(x, cqc_data){
  invisible(apply(x, 1, function(row){
    sn <- row[["FCS"]]
    cf <- cqc_data[[sn]]
    cf_rename_marker(cf, row[["from"]], row[["to"]])
  })
  )
}

#' @export
cqc_remove_not_in_reference <- function(x, ...)UseMethod("cqc_remove_not_in_reference")
#' @export
cqc_remove_not_in_reference.cqc_report_channel <- function(x, cqc_data){
  for(sn in names(x))
  {
    check_result <- x[[sn]]
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]
    if(length(unknown) > 0 && length(missing) == 0)
    {
      cf <- cqc_data[[sn]]
      j <- match(unknown, colnames(cf))
      subset_cytoframe_by_cols(cf@pointer, j - 1)

    }
  }
}
#' @export
cqc_remove_not_in_reference.cqc_report_marker <- function(x, cqc_data){
  for(sn in names(x))
  {
    check_result <- x[[sn]]
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]
    if(length(unknown) > 0 && length(missing) == 0)
    {
      cf <- cqc_data[[sn]]
      cf_rename_marker(cf, unknown, "")

    }
  }
}
#' @export
cqc_drop_samples <- function(cqc_data, check_results){

    cqc_data[-match(names(check_results), names(cqc_data))]
}

#' @export
#' @importFrom dplyr filter arrange pull mutate group_indices
#' @importFrom tidyr separate
cqc_group_by_panel <- function(cqc_data, delimiter = "|"){
  sep <- paste0(delimiter, delimiter)#double delimiter for sep params and single delimiter for sep channel and marker
  keys <- sapply(cqc_data, function(cf){
              cf_get_params_tbl(cf) %>%
                arrange(channel) %>%
                  filter(is.na(marker) == FALSE) %>%
                    unite(panel, channel, marker, sep = delimiter) %>%
                      pull(panel) %>%
                        paste(collapse = sep)
  })
  res <- tibble(FCS = names(keys), panel = keys)

  res <- res %>% mutate(panel_id = group_indices(res, panel))

  #
  # res <- strsplit(res, split= sep, fixed = "TRUE")[[1]]
  # res <- tibble(reference = res)
  # if(type == "marker")
  #   res <- separate(res, channel, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
  class(res) <- c("cqc_group", class(res))
  class(res) <- c("cqc_group_panel", class(res))
  attr(res, "delimiter") <- delimiter
  res
}

#' @export
summary.cqc_group_panel <- function(object){
  object %>% group_by(panel)
}

#' @importFrom tidyr separate_rows
#' @importFrom dplyr distinct count
cqc_get_panel_info <- function(x){
  delimiter <- attr(x, "delimiter")
  sep <- paste0(delimiter, delimiter)
  x %>% count(panel_id, panel) %>%
      rename(nFCS = n) %>%
      separate_rows(panel, sep = paste0("\\Q", sep, "\\E")) %>%
        separate(panel, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
}

#' @export
format.cqc_group_panel <- function(x){

  cqc_get_panel_info(x) %>%
    kable() %>%
     kable_styling("bordered", full_width = F, position = "left") %>%
      collapse_rows(columns = c(1,4), "top")%>%
        row_spec(0, background = "#e5f5e0", color = "black")
}

#' @importFrom dplyr group_split
#' @export
diff.cqc_group_panel <- function(x){
  cqc_get_panel_info(x) %>%
    select(panel_id, channel, marker) %>%
      group_split(panel_id)
}


#' @export
split.cqc_group_panel <- function(x, cqc_data){

}
