
#' @export
cqc_fix <- function(x, ...)UseMethod("cqc_fix")

#' Apply the cqc_solution
#'
#' Peform the actual fixing action (i.e update or delete)
#' @param x the cqc_solution returned by 'find_solution' calls
#' @param ... addiitional arguments not for the user.
#' @importFrom dplyr rowwise do
#' @export
cqc_fix.cqc_solution <- function(x,...){
  group <- attr(x, "group")
  cqc_data <- attr(group, "data")
  type <- sub("cqc_solution_", "", class(x)[[1]])
  invisible(group %>% inner_join(x, "group_id") %>%
              select(object, from, to) %>% distinct() %>% rowwise() %>% do({
                obj <- cqc_data[[.[["object"]]]]
                if(is.na(.[["to"]]))
                {
                  cqc_delete(obj, .[["from"]], type)

                }else
                  cqc_update(obj, .[["from"]], .[["to"]], type)
                tibble()
              }))

}
#' @export
cqc_delete <- function(x, ...)UseMethod("cqc_delete")
#' @export
cqc_delete.cytoframe <- function(x, value, type, ...){
  if(type == "channel")
  {
    cols <- flowWorkspace::colnames(x)
    j <- which(!cols %in% value)
    flowWorkspace:::subset_cytoframe_by_cols(x@pointer, j - 1)
  }else if(type == "marker")
  {
    cf_rename_marker(x, value, "")
  }else if(type == "keyword")
  {
    cf_keyword_delete(x, value)

  }else
    stop("don't know how to proceed!")
}
#' @export
cqc_delete.GatingSet <- function(x, ...){
  cs <- gs_cyto_data(x)#cs is a new view
  lapply(cs, cqc_delete, ...) # deleting channel is done through subsetting, thus only affect the view not the original data
  gs_cyto_data(x) <- cs# thus need to assign it back to gs to take effect

}
#' @export
cqc_update <- function(x, ...)UseMethod("cqc_update")
#' @export
cqc_update.cytoframe <- function(x, from, to, type,...){
  if(type == "channel")
  {
    flowWorkspace:::setChannel(x@pointer, from, to)
  }else if(type == "marker")
  {
    cf_rename_marker(x, from, to)
  }else if(type == "keyword")
  {
    cf_keyword_rename(x, from, to)

  }else
    stop("don't know how to proceed!")

}
#' @export
cqc_update.GatingSet <- function(x, from, to, type,...){
  if(type == "channel")
    gs_update_channels(x, map = data.frame(old = from
                                           , new = to
                                            )
                        )
  else
  {
   cs <-  gs_cyto_data(x)
   lapply(cs, cqc_update, from, to, type)#cs point to the original data, no need to assigning it back

  }
}
