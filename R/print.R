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

  x <- x %>%
    mutate("Proposed change" = ifelse(is.na(to), cell_spec(from, strikeout = TRUE), paste(from, to, sep = " --> "))) %>%
    select(-c(from, to)) %>%
    kable(escape = T) %>%
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
  type <- sub("cqc_group_", "", class_names[2])
  if (type != "panel") {
    type <- as.symbol(type)
    x <- group_by(x, group_id, nObject) %>%
      summarise(!!type := paste(!!type, collapse = ", ")) %>%
      arrange(desc(nObject))
    class(x) <- class_names
  }


  x
}
#' #' @export
#' print.cqc_group_summary <- function(x, collapse = TRUE, ...){
#'   if(collapse)
#'     x <- collapse_params(x)
#'   class(x) <- class(x)[-1]
#'   print(x)
#' }

#' @importFrom dplyr ungroup everything
#' @export
knit_print.cqc_group_summary <- function(x, collapse = TRUE, ...) {
  n <- nrow(x)
  if (collapse) {
    x <- collapse_params(x)
  }

  collaspse_idx <- match("group_id", colnames(x))
  if (is(x, "cqc_group_panel")) {
    collaspse_idx <- c(collaspse_idx, match("nObject", colnames(x)))
  }

  x <- x %>%
    arrange(desc(nObject)) %>%
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
