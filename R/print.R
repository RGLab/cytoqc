#' @export
print.cqc_cf_list <- function(x, ...) {
  cat("cytoqc data: \n")
  cat(length(x), " samples \n")
}

#' @importFrom knitr knit_print
#' @export
knit_print.cqc_reference <- function(x, ...) {
  x %>%
    kable() %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 12) %>%
    row_spec(0, background = "gray", color = "black") %>%
    knit_print()
}

#' @export
print.cqc_match_result <- function(x, ...) {
  print(as_tibble(x))
}

#' @export
print.cqc_match_result_and_solution <- function(x, ...) {
  print(format(x, ...))
}


#' @export
knit_print.cqc_match_result_and_solution <- function(data, ...) {
  datatable(format(data, ...))
}
#' @export
#' @importFrom purrr map_dfc
format.cqc_match_result_and_solution <- function(x, ...) {
  #combine match res and recommened solution to present it as wide format for easy viewing the data
  ref <- x[["ref"]]
  match_result <- x[["match_result"]]
  gids <- names(x[["match_result"]])
  tbl <- map_dfc(gids, function(gid) {
      df <- filter(x[["solution"]], group_id == gid)
      this_res <- match_result[[gid]]
      #init column
      col_to_show <- rep("*", length(ref))

      if(!is.null(this_res))
      {
        missing <- this_res[["missing"]]
        unknown <- this_res[["unknown"]]
         #drop the recommended deletion
        df1 <- filter(df, !is.na(to))
        #fill the refs with the recommended edit
        matched.ref <- df1[["to"]]
        matched.target <- df1[["from"]]
        idx <- match(matched.ref, ref)
        col_to_show[idx] <- matched.target
        # browser()
        #fill the unmatched refs
        unmatched.ref <- missing[!missing%in%matched.ref]
        col_to_show[match(unmatched.ref, ref)] <- NA
        #append the unmatched target
        if(length(matched.target)>0)
          unmatched.target <- unknown[-match(matched.target, unknown)]#exclude the matched ones
        else
          unmatched.target <- unknown
        if(length(unmatched.target) == 0)
          unmatched.target <- ""
        else
          unmatched.target <- paste(unmatched.target, collapse = ",")
        col_to_show <- c(col_to_show, unmatched.target)

        ##append the redundant item
        torm <- filter(df, is.na(to))[["from"]]
        if(length(torm) == 0)
          torm <- ""
        else
          torm <- paste(torm, collapse = ",")
        col_to_show <- c(col_to_show, torm)
      }else
        col_to_show <- c(col_to_show, "", "")

      col_to_show <- list(col_to_show)
      names(col_to_show) <- gid
      col_to_show
    })
  tbl <- cbind(c(ref, "", ""), tbl)
  tbl <- cbind(c(rep("", length(ref)), "Unmatched", "To Delete"), tbl)
  colnames(tbl)[1:2] <- c("", "Ref")
  #rm last two rows if they are all empty
  torm <- tbl[nrow(tbl),-1]
  if(isTRUE(all(torm == "")))
    tbl <- tbl[-nrow(tbl),]
  unmatched <- tbl[nrow(tbl), -1]
  if(isTRUE(all(unmatched == "")))
    tbl <- tbl[-nrow(tbl),]
  tbl
}
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling row_spec
#' @export
knit_print.cqc_match_result <- function(x, ...) {
  if (length(x) == 0) {
    res <- kable(data.frame(object = "All passed"), col.names = NULL) %>% row_spec(1, color = "green")
  } else {
    res <- kable(as_tibble(x)) %>%
      row_spec(0, background = "#9ebcda", color = "black")
  }

  res %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 12) %>%
    column_spec(1, bold = TRUE) %>%
    knit_print()
}

#' @importFrom htmltools htmlEscape
#' @importFrom kableExtra collapse_rows cell_spec column_spec
#' @importFrom dplyr %>% select distinct mutate_if
#' @importFrom tidyr unite
#' @export
knit_print.cqc_solution <- function(x, itemize = FALSE, ...) {
  if (!itemize) {
    x <- x %>%
      select(-1) %>%
      distinct()
  }

  x <- x %>% mutate(from = htmlEscape(from), to = htmlEscape(to)) %>%
    mutate("Proposed change" = ifelse(is.na(to), cell_spec(from, strikeout = TRUE), paste(from, to, sep = " --> "))) %>%
    select(-c(from, to)) %>%
    kable(escape = F) %>%
    kable_styling(c("bordered", "condensed"), full_width = F, position = "left", font_size = 12) %>%
    collapse_rows(columns = 1, "top") %>%
    row_spec(0, background = "#e5f5e0", color = "black")
  if (itemize) {
    x <- x %>% column_spec(1, bold = TRUE)
  }
  knit_print(x)
}
#' @importFrom dplyr summarise
collapse_params <- function(x, ...) {
  class_names <- class(x)
  type <- sub("cqc_check_", "", class_names[2])
  if (type != "panel") {
    type <- as.symbol(type)
    x <- group_by(x, group_id, nObject) %>%
      summarise(!!type := paste(!!type, collapse = ", ")) %>%
      arrange(desc(nObject))
    class(x) <- class_names
  }


  x
}

#' @export
print.cqc_check <- function(x, ...){
   print(summary(x), ...)
}

object_type <- function(x){
  dat <- attr(x, "data")
  if(is(dat, "cqc_gs_list"))
    "GatingSet"
  else
    "FCS"
}

#' @export
print.cqc_check_summary <- function(x, collapse = TRUE, ...){
  if(collapse)
    x <- collapse_params(x)
  type <- object_type(x)
  colnames(x)[match("nObject", colnames(x))] = paste0("n", type)
  class(x) <- class(x)[-(1:3)]
  print(x)
}

#' @export
knit_print.cqc_check <- function(x, ...){
  knit_print(summary(x), ...)
}
#' Customized knit print for cqc_check_summary
#'
#'
#' @param x cqc_check_summary object returned by 'summary' call on `cqc_check`
#' @param collapse whether to collapse the same information within each group
#' @param ... not used
#' @importFrom dplyr ungroup everything
#' @export
knit_print.cqc_check_summary <- function(x, collapse = TRUE, ...) {
  n <- nrow(x)
  if (collapse) {
    x <- collapse_params(x)
  }
  type <- object_type(x)
  type <-  paste0("n", type)
  colnames(x)[match("nObject", colnames(x))] = type

  collaspse_idx <- match("group_id", colnames(x))
  if (is(x, "cqc_check_panel")) {
    collaspse_idx <- c(collaspse_idx, match(type, colnames(x)))
  }

  x <- x %>%
    arrange(desc(get(type))) %>%
    kable() %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 12)
  # browser()
  if (n > 0) {
    x <- x %>%
      collapse_rows(columns = collaspse_idx, "top") %>%
      row_spec(0, background = "#e5f5e0", color = "black")
  }

  knit_print(x)
}

#' @export
#' @importFrom tidyr spread
format.cqc_check_panel <- function(x, anchor = c("channel", "marker"), ...){
  # x <- summary(x)
  anchor <- match.arg(anchor)
  if(anchor == "channel")
    value <- "marker"
  else
    value <- "channel"
  #long to wide
  x %>% summary %>%
    mutate(group_id := paste("group", group_id), nObject := paste0("(n=", nObject, ")")) %>%
    unite(grp, group_id, nObject, sep = "") %>% #merge grp cols
    spread(grp, !!value) %>%
    `class<-`(value = c("cqc_check_panel_wide", class(x)))
}

#' @export
print.cqc_check_panel_wide <- function(x, ...){
    x %>% `class<-`(value = class(x)[-(1:3)]) %>%
    print
}
#' @export
print.cqc_check_panel <- function(x, ...){
  format(x, ...) %>%
    print
}

#' @export
knit_print.cqc_check_panel <- function(x, ...){

  format(x, ...) %>%
    knit_print
}
#' @export
#' @importFrom DT datatable
knit_print.cqc_check_panel_wide <- function(x, ...){
  x%>%
    `class<-`(value = class(x)[-(1:3)])  %>% datatable %>% knit_print
}


