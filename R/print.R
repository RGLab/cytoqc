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
  kable_styling("bordered", full_width = F, position = "left") %>%
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
    kable_styling("bordered", full_width = F, position = "left") %>%
    column_spec(1, bold = TRUE) %>%
    knit_print

}

#' @importFrom kableExtra collapse_rows
#' @importFrom dplyr %>% select distinct
#' @importFrom tidyr unite
#' @export
knit_print.cqc_solution <- function(x, itemize = FALSE, ...){
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
  knit_print(x)

}

#' @export
knit_print.cqc_group_panel_summary <- function(x, ...){

  x %>%
    kable() %>%
     kable_styling("bordered", full_width = F, position = "left") %>%
      collapse_rows(columns = c(1,4), "top")%>%
        row_spec(0, background = "#e5f5e0", color = "black") %>%
    knit_print

}
