#' @export
print.cqc_data <- function(x){
  cat("cytoqc data: \n")
  cat(length(x), " samples \n")
}

#' @importFrom knitr knit_print
#' @export
knit_print.cqc_reference <- function(x, ...){
  x %>%
  kable() %>%
  kable_styling("bordered", full_width = F, position = "left", font_size = 10) %>%
  row_spec(0, background = "gray", color = "black") %>%
  knit_print

}

#' @export
print.cqc_report_params <- function(x, ...){
  print(as_tibble(x))
}


#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @export
knit_print.cqc_report_params <- function(x, ...){
  if(length(x) == 0)
    res <- kable(data.frame(FCS = "All passed"), col.names = NULL)%>% row_spec(1, color = "green")
  else
    res <- kable(as_tibble(x)) %>%
      row_spec(0, background = "#9ebcda", color = "black")

  res %>%
    kable_styling("bordered", full_width = F, position = "left", font_size = 10) %>%
    column_spec(1, bold = TRUE) %>%
    knit_print

}

#' @importFrom kableExtra collapse_rows cell_spec
#' @importFrom dplyr %>% select distinct mutate_if
#' @importFrom tidyr unite
#' @export
knit_print.cqc_solution <- function(x, itemize = FALSE, ...){
  if(!itemize)
    x <- x %>% select(-1) %>% distinct()

  x <- x %>%
    mutate("Proposed change" = ifelse(is.na(to), cell_spec(from, strikeout = TRUE), paste(from, to, sep = " --> "))) %>%
    select(-c(from, to)) %>%
    kable() %>%
    kable_styling(c("bordered", "condensed"), full_width = F, position = "left", font_size = 10) %>%
    collapse_rows(columns = 1, "top")%>%
    row_spec(0, background = "#e5f5e0", color = "black")
  if(itemize)
    x <- x %>% column_spec(1, bold = TRUE)
  knit_print(x)

}
#' @export
print.cqc_group_summary <- function(x, collapse = TRUE, ...){
  collapse_params(x)
}
#' @export
knit_print.cqc_group_summary <- function(x, collapse = TRUE, ...){
  cn <- colnames(x)
  idx <- cn %in% c("group_id", "nFCS")
  #reorder column to place groupid first and nFcs last
  x %>% select(c(which(idx)[1], which(!idx), which(idx)[2])) %>%
    arrange(desc(nFCS)) %>%
    kable() %>%
     kable_styling("bordered", full_width = F, position = "left", font_size = 10) %>%
      collapse_rows(columns = c(1, length(idx)), "top")%>%
        row_spec(0, background = "#e5f5e0", color = "black") %>%
    knit_print

}
