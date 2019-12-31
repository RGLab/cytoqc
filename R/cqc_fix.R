
#' @export
cqc_fix <- function(x, ...)UseMethod("cqc_fix")

#' Apply the cqc_solution
#'
#' Peform the actual fixing action (i.e update or delete)
#' @param x the cqc_solution returned by 'find_solution' calls
#'
#' @importFrom dplyr rowwise do
cqc_fix.cqc_solution <- function(x, func){
  group <- attr(x, "group")
  cqc_data <- attr(group, "data")

  invisible(group %>% inner_join(x, "group_id") %>%
              select(FCS, from, to) %>% distinct() %>% rowwise() %>% do({
                cf <- cqc_data[[.[["FCS"]]]]
                if(is.na(.[["to"]]))
                {
                  if(is(x, "cqc_solution_channel"))
                  {
                    cols <- flowWorkspace::colnames(cf)
                    j <- which(!cols %in% .[["from"]])
                    flowWorkspace:::subset_cytoframe_by_cols(cf@pointer, j - 1)
                  }else if(is(x, "cqc_solution_marker"))
                  {
                    cf_rename_marker(cf, .[["from"]], "")
                  }else
                    stop("don't know how to proceed!")

                }else
                  func(cf, .[["from"]], .[["to"]])
                tibble()
              }))

}
#' @export
cqc_fix.cqc_solution_channel <- function(x){
  cqc_fix.cqc_solution(x, function(cf,...)flowWorkspace:::setChannel(cf@pointer, ...))
}

#' @export
cqc_fix.cqc_solution_marker <- function(x){
  cqc_fix.cqc_solution(x, cf_rename_marker)

}

#' Helper function to remove the outlier groups that can't be fixed
#'
#' @param cqc_groups the object returned by 'cqc_check'
#' @param id the group id to be dropped from the dataset
#' @export
cqc_drop_groups <- function(cqc_groups, id){
  cqc_data <- attr(cqc_groups, "data")
  torm <- filter(cqc_groups, group_id == id) %>% pull(FCS) %>% unique()
  cqc_data <- cqc_data[-match(torm, names(cqc_data))]
  class(cqc_data) <- "cqc_data"
  cqc_data
}



# cqc_get_panel_info <- function(x){
#   delimiter <- attr(x, "delimiter")
#   sep <- paste0(delimiter, delimiter)
#   x %>% count(panel_id, panel) %>%
#       rename(nFCS = n) %>%
#       separate_rows(panel, sep = paste0("\\Q", sep, "\\E")) %>%
#         separate(panel, c("channel", "marker"), sep = paste0("\\Q", delimiter, "\\E"))
# }

