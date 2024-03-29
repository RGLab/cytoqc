#' Apply the cqc_solution
#'
#' Peform the actual fixing action (i.e update or delete)
#' @param x the \code{cqc_solution} returned by \code{\link{cqc_match}} calls
#' @param ... addiitional arguments not for the user.
#' @examples
#' # Read in FCS files with inconsistencies
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#'
#' # Check for marker inconsitencies
#' groups <- cqc_check(qc_cf_list, type = "marker")
#'
#' # Attempt to fix them automatically
#' match_result <- cqc_match(groups, ref = c("CD14 PerCP", "CD15 FITC", "CD33 APC", "CD45 PE", "FSC-Height", "SSC-Height", "Time"))
#'
#' # Add a manual match that automatic matching could not find
#' match_result <- cqc_match_update(match_result, map = c("PTPRC PE" = "CD45 PE"))
#'
#' # Apply the fix to the original cytoframes
#' cqc_fix(match_result)
#' @export
cqc_fix <- function(x, ...) UseMethod("cqc_fix")

#' @export
cqc_fix.default <- function(x, ...) {
  stop("The input is not a valid 'cqc_match' result!\nPlease make sure to follow the right order of the 'cqc_check-->cqc_match-->cqc_fix' workflow!")
}

#' @export
cqc_fix.cqc_match_result_and_solution <- function(x, ...) {
  cqc_fix(x[["solution"]])
}
#' @importFrom dplyr rowwise do
#' @export
cqc_fix.cqc_solution <- function(x, ...) {
  group <- attr(x, "group")
  cqc_data <- attr(group, "data")
  type <- sub("cqc_solution_", "", class(x)[[1]])
  invisible(group %>% inner_join(x, "group_id") %>%
    select(object, from, to) %>% distinct() %>% rowwise() %>% do({
      obj <- cqc_data[[.[["object"]]]]
      if (!is.na(.[["from"]])&&is.na(.[["to"]])) {#xxx -> NA
        cqc_delete(obj, .[["from"]], type)
      } else if(!is.na(.[["from"]])&&!is.na(.[["to"]])){#xxx -> yyy
        cqc_update(obj, .[["from"]], .[["to"]], type)
      } else if(is.na(.[["from"]])&&!is.na(.[["to"]]))#NA -> yyy
      {
        cqc_insert(obj, .[["to"]], type)
      }else#NA -> NA
        stop("both 'from' and 'to' are NA!")
      tibble()
    }))
}

#' insertion methods for cyto data
#'
#' It is typically called automatically by cqc_fix call
#'
#' @param x cytoframe or GatingSet
#' @param value the value to be inserted
#' @param type one of qc task only  "keyword" is valid for now.
#' @param ... unused
cqc_insert <- function(x, value, type, ...) UseMethod("cqc_insert")

cqc_insert.cytoframe <- function(x, value, type, ...) {
  if  (type == "keyword") {
    cf_keyword_insert(x, value, "")
  } else {
    stop("insert ", type, " not supported yet!")
  }
}

cqc_insert.GatingSet <- function(x, value, type, ...) {
  cs <- gs_cyto_data(x) # cs is a new view
  lapply(cs, function(cf) {
    cqc_insert(cf, value, type, ...)
  })

}


#' Delete methods for cyto data
#'
#' It is typically called automatically by cqc_fix call
#'
#' @param x cytoframe or GatingSet
#' @param value the value to be deleted
#' @param type one of qc task "channel", "marker", "keyword", "gate"
#' @param ... unused
#' @export
cqc_delete <- function(x, value, type, ...) UseMethod("cqc_delete")
#' @export
cqc_delete.cytoframe <- function(x, value, type, ...) {
  if (type == "channel") {
    cols <- flowWorkspace::colnames(x)
    j <- which(!cols %in% value)
    flowWorkspace:::subset_cytoframe_by_cols(x@pointer, as.integer(j - 1))
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

#' Update methods for cyto data
#'
#' It is typically called automatically by cqc_fix call
#'
#' @param x cytoframe or GatingSet
#' @param from the old value to be updated
#' @param to the old value to be updated. Alternatively when type is "marker" it can be a named character which provides
#'          a pair of channel vs marker and 'from' is ignored , instead 'channel' is used as the reference to select target entry to update
#'
#' @param type one of qc task "channel", "marker", "keyword", "gate"
#' @param ... unused
#' @export
cqc_update <- function(x, from, to, type, ...) UseMethod("cqc_update")
#' @export
cqc_update.cytoframe <- function(x, from, to, type, ...) {
  if (type == "channel") {
    flowWorkspace:::setChannel(x@pointer, from, to)
  } else if (type == "marker") {
    if(length(names(to)) == 0)
    {
      if(from == "")
        stop("Can't rename the original empty marker directly!Please specify 'to' as a named character which provides
                  a pair of channel vs marker to select target marker to update")
      cf_rename_marker(x, from, to)
    }else
    {
      if(from != "")
        warning("'from' argument is ignored since 'to' is a named character which provides
                  a pair of channel vs marker to select target marker to update")
      markernames(x) <- to
    }

  } else if (type == "keyword") {
    cf_keyword_rename(x, from, to)
  } else {
    stop("don't know how to proceed!")
  }
}
#' @export
cqc_update.GatingSet <- function(x, from, to, type, ...) {
  is_gh <- is(x, "GatingHierarchy")
  if (type == "channel") {
    gs_update_channels(x, map = data.frame(
      old = from,
      new = to
    ))
  } else if (type == "gate") {
    if(is_gh)
      gh_pop_set_name(x, from, to)
    else
      gs_pop_set_name(x, from, to)
  } else {
    if(is_gh){
      cf <- gh_pop_get_data(x, returnType = "cytoframe")
      cqc_update(cf, from, to, type)
    }else{
      cs <- gs_cyto_data(x)
      lapply(cs, cqc_update, from, to, type) # cs point to the original data, no need to assigning it back
    }
  }
}
#' @export
cqc_fix.cqc_match_result_panel <- function(x, ...) {
  ref <- attr(x, "ref")
  by <- attr(x, "by")
  cqc_fix_panel(x, ref , by)
}

#' Update panel info
#'
#' It updates either 'channel' or 'marker' info depends on the value of 'ref.col'
#' if ref.col is set to 'channel', then 'channel' column from 'panel' is used as the
#' reference to be matched against 'cqc_data' and 'marker'  will be set according to 'panel'
#'
#' @param x "cqc_check_panel" result from 'cqc_check_panel' or 'cqc_check'call
#' @param ref the reference group id or a custom panel that is represented as tibble contains channel and marker columns (typically returned by 'cf_get_panel` )
#' @param by the column used as the fixed anchor, either 'channel' or 'marker'
#' @param select the groups to be updated. Default is all
#' @examples
#' \dontrun{
#'      res <- cqc_check(cqc_data, "panel")
#'      cqc_fix_panel(res, 1, "channel") #get panel info from group 1 and set markers of other groups  by using channel as reference
#'
#' }
cqc_fix_panel <- function(x, ref, by, select = NULL)
{
  stopifnot(is(x, "cqc_check_panel"))
  data <- attr(x, "data")

  if(is(ref, "numeric"))
  {
    #extract ref panel
    pnl <- summary(x) %>%
            filter(group_id == ref) %>%
            select(channel, marker)
    #exclude ref group from targets
    x <- filter(x,group_id != ref)
    ref <- pnl
  }
  if (!is.null(select)) {
    x <- filter(x, group_id %in% select)#select the target groups
  }
  obj <- unique(x[["object"]])
  for(i in data[obj])
    cqc_set_panel(i, ref, by)
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
#' @param ... unused
cqc_set_panel <- function(x, panel, ref.col, ...) UseMethod("cqc_set_panel")

cqc_set_panel.cqc_data <- function(x, ...){
  for(i in x)
    cqc_set_panel(i, ...)
}

cqc_set_panel.GatingSet <- function(x,  panel, ref.col, ...){
  cols <- c("channel", "marker")
  ref.col <- match.arg(ref.col, cols)
  cs <- gs_cyto_data(x)
  is_gh <- is(x, "GatingHierarchy")
  if(is_gh)
    cf <- gh_pop_get_data(x, returnType = "cytoframe")
  else
    cf <- get_cytoframe_from_cs(cs, 1)

  if (ref.col == "marker") {
    if(!setequal(colnames(panel), cols))
      stop("invalid 'panel' info")

    target.col <- cols[-match(ref.col, cols)]
    old <- paste0(target.col, ".x")
    new <- paste0(target.col, ".y")

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
    if(is_gh)
      cqc_set_panel(cf, panel, ref.col, ...)
    else
      lapply(cs, cqc_set_panel,  panel, ref.col, ...) # cs point to the original data, no need to assigning it back
  }
}
cqc_set_panel.cytoframe <- function(x, panel, ref.col, ...){
  cols <- c("channel", "marker")
  ref.col <- match.arg(ref.col,cols)
  if(!setequal(colnames(panel), cols))
    stop("invalid 'panel' info")
  target.col <- cols[-match(ref.col, cols)]
  old <- paste0(target.col, ".x")
  new <- paste0(target.col, ".y")

  ref.markers <- panel[["marker"]]
  ref.channels <- panel[["channel"]]
  pnl1 <- cf_get_panel(x, skip_na = FALSE)
  if(ref.col == 'channel')
  {
    names(ref.markers) <- ref.channels
    markernames(x) <- ref.markers
  }else
  {
    markers <- pnl1[["marker"]]

    #find idx of ref in target
    idx <- match(ref.markers, markers)
    colnames(x)[idx] <- ref.channels
  }
}
