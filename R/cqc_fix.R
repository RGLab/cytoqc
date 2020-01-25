
#' @export
cqc_fix <- function(x, ...) UseMethod("cqc_fix")

#' Apply the cqc_solution
#'
#' Peform the actual fixing action (i.e update or delete)
#' @param x the cqc_solution returned by 'find_solution' calls
#' @param ... addiitional arguments not for the user.
#' @importFrom dplyr rowwise do
#' @export
cqc_fix.cqc_solution <- function(x, ...) {
  group <- attr(x, "group")
  cqc_data <- attr(group, "data")
  type <- sub("cqc_solution_", "", class(x)[[1]])
  invisible(group %>% inner_join(x, "group_id") %>%
    select(object, from, to) %>% distinct() %>% rowwise() %>% do({
      obj <- cqc_data[[.[["object"]]]]
      if (is.na(.[["to"]])) {
        cqc_delete(obj, .[["from"]], type)
      } else {
        cqc_update(obj, .[["from"]], .[["to"]], type)
      }
      tibble()
    }))
}
#' @export
cqc_delete <- function(x, ...) UseMethod("cqc_delete")
#' @export
cqc_delete.cytoframe <- function(x, value, type, ...) {
  if (type == "channel") {
    cols <- flowWorkspace::colnames(x)
    j <- which(!cols %in% value)
    flowWorkspace:::subset_cytoframe_by_cols(x@pointer, j - 1)
  } else if (type == "marker") {
    cf_rename_marker(x, value, "")
  } else if (type == "keyword") {
    cf_keyword_delete(x, value)
  } else {
    stop("don't know how to proceed!")
  }
}
#' @export
cqc_delete.GatingSet <- function(x, value, type, ...) {
  cs <- gs_cyto_data(x) # cs is a new view
  if (type == "channel") {
    # deleting channel is done through subsetting, thus only affect the view not the original data
    #    get_cytoframe_from_cs is getting a new view  of cf thus can't be updated inplace
    cols <- flowWorkspace::colnames(x)
    j <- which(!cols %in% value)
    cs <- cs[, j]
    gs_cyto_data(x) <- cs # thus need to assign it back to gs to take effect
  } else if (type == "gate") {
    gs_pop_set_visibility(x, value, FALSE)
  } else {
    lapply(cs, function(cf) {
      cqc_delete(cf, value, type)
    })
  }
}
#' @export
cqc_update <- function(x, ...) UseMethod("cqc_update")
#' @export
cqc_update.cytoframe <- function(x, from, to, type, ...) {
  if (type == "channel") {
    flowWorkspace:::setChannel(x@pointer, from, to)
  } else if (type == "marker") {
    cf_rename_marker(x, from, to)
  } else if (type == "keyword") {
    cf_keyword_rename(x, from, to)
  } else {
    stop("don't know how to proceed!")
  }
}
#' @export
cqc_update.GatingSet <- function(x, from, to, type, ...) {
  if (type == "channel") {
    gs_update_channels(x, map = data.frame(
      old = from,
      new = to
    ))
  } else if (type == "gate") {
    gs_pop_set_name(x, from, to)
  } else {
    cs <- gs_cyto_data(x)
    lapply(cs, cqc_update, from, to, type) # cs point to the original data, no need to assigning it back
  }
}

#' Update panel info
#'
#' It updates either 'channel' or 'marker' info depends on the value of 'ref.col'
#' if ref.col is set to 'channel', then 'channel' column from 'panel' is used as the
#' reference to be matched against 'cqc_data' and 'marker'  will be set according to 'panel'
#'
#' @param x cqc_data
#' @param panel a tibble contains channel and marker colums (typically returned by 'cf_get_panel` )
#' @param ref.col the column used as the reference, either 'channel' or 'marker'
#' @export
cqc_set_panel <- function(x, panel, ref.col, ...) UseMethod("cqc_set_panel")

#' @export
cqc_set_panel.cqc_data <- function(x, ...){
  for(i in x)
    cqc_set_panel(i, ...)
}

#' @export
cqc_set_panel.GatingSet <- function(x,  panel, ref.col, ...){
  cols <- c("channel", "marker")
  ref.col <- match.arg(ref.col, cols)
  cs <- gs_cyto_data(x)

  if (ref.col == "marker") {
    if(!setequal(colnames(panel), cols))
      stop("invalid 'panel' info")

    target.col <- cols[-match(ref.col, cols)]
    old <- paste0(target.col, ".x")
    new <- paste0(target.col, ".y")
    cf <- get_cytoframe_from_cs(cs, 1)

    tbl <- cf_get_panel(cf, skip_na = TRUE) %>%
      inner_join(panel, by = ref.col) %>%
      filter(get(old) != get(new))
    if(nrow(tbl) > 0)
      gs_update_channels(x, map = data.frame(
                                        old = tbl[[old]],
                                        new = tbl[[new]]
                                      )
                       )
  } else {
    lapply(cs, cqc_set_panel,  panel, ref.col, ...) # cs point to the original data, no need to assigning it back
  }
}
#' @export
cqc_set_panel.cytoframe <- function(x, panel, ref.col, ...){
  cols <- c("channel", "marker")
  ref.col <- match.arg(ref.col,cols)
  if(!setequal(colnames(panel), cols))
    stop("invalid 'panel' info")
  target.col <- cols[-match(ref.col, cols)]
  old <- paste0(target.col, ".x")
  new <- paste0(target.col, ".y")
  # fun <- get(paste0("cf_rename_", target.col))
  cf_get_panel(x, skip_na = TRUE) %>%
    inner_join(panel, by = ref.col) %>%
    filter(get(old) != get(new)) %>%
    rowwise() %>% do({
      # browser()
      # fun(x, .[[old]], .[[new]])
      cqc_update(x, .[[old]], .[[new]], type = target.col)
    })

}
