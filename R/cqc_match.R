#' find the the difference between the reference and target group
#'
#'
#' @param x \code{cqc_check} result returned by cqc_check call
#' @param ...
#'
#'        ref -- specifies the reference, which can be either an integer group id or a character vector giving the actual values of the reference
#'
#'        select -- the group ids selected for processing
#'
#'        type -- the qc type (either "channel", "marker", "gate"), automatically determined by the type of \code{x}
#'
#'        max.distance -- Maximum distance allowed for a match. This is passed to the max.distance argument in \code{\link{agrep}}.
#'
#'        partial whether -- to do the partial sub string matching before the approximate string matching
#'
#' @examples
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' channel_groups <- cqc_check(qc_cf_list, type = "channel")
#' summary(channel_groups)
#'
#' # Match to a reference group of samples
#' channel_match <- cqc_match(channel_groups, 3)
#'
#' # Match to a character vector of target channel names
#' channel_match <- cqc_match(channel_groups, ref = c("FL1-H", "FL2-A", "FL2-H", "FL3-H", "FL4-H", "FSC-Height", "SSC-Height", "Time"))
#' @export
cqc_match <- function(x, ...) UseMethod("cqc_match")

#' @export
cqc_match.default <- function(x, ...) {
  stop("The input is not a valid 'cqc_check' result!\nPlease make sure to follow the right order of 'cqc_check-->cqc_match-->cqc_fix' workflow!")
}

#' @export
cqc_match.cqc_check_channel <- function(x, ...) {
  res <- match_reference(x, type = "channel", ...)
  res
}

#' @export
cqc_match.cqc_check_marker <- function(x, ...) {
  res <- match_reference(x, type = "marker", ...)
  res
}

#' @export
cqc_match.cqc_check_keyword <- function(x, ...) {
  res <- match_reference(x, type = "keyword", ...)
  res
}
#' @export
cqc_match.cqc_check_gate <- function(x, ...) {
  res <- match_reference(x, type = "gate", ...)
  res
}

#' @export
#' @importFrom dplyr group_walk
cqc_match.cqc_check_panel <- function(x, ref, ...) {
  by <- attr(x, "by")
  #check if anchor is already standardized
  stopifnot(is(ref, "numeric"))
  # Drop rows with NA in anchor
  x <- x %>% filter(!is.na(!!as.name(by)))
  
  ref_by <- filter(x, group_id == ref)[[by]]
  x %>% filter(group_id != ref) %>%
    group_by(group_id) %>% group_walk(function(df,...){
      if(!setequal(df[[by]], ref_by))
        stop(by, " is not consistent across panel groups!Please standardize it first!")
    })
  #simply store the ref and by(or anchor)
  #the actual matching(or alignment) is done in format method
  attr(x, "ref") <- ref
  attr(x, "by") <- by
  class(x) <- c("cqc_match_result_panel", "cqc_match_result", class(x))
  x
}
summary.cqc_match_result_panel <- function(object, ...) {
  cls <- class(object)[1:2]
  class(object) <- class(object)[-c(1:2)]
  x <- summary.cqc_check(object, ...)
  x
}
#' find the the difference between the reference and target group
#'
#' Only used internally.
#'
#' @param x cqc report returned by set_reference call
#' @param ref specifies the reference, which can be either an integer group id or a characte vector giving the actual values of the reference
#' @param select the group ids selected for processing
#' @param type the qc type (either "channle", "marker")
#' @importFrom dplyr bind_rows group_keys group_by
#' @importFrom purrr set_names
#' @noRd
match_reference <- function(x, ref, select = NULL, type, delimiter = "|", ...) {
  res <- summary(x)#get group-wise check report
  if (is.numeric(ref)) {#fetch reference vector from check result if ref is an id
    refid <- ref
    ref <- res %>%
      filter(group_id %in% ref) %>%
      pull(type)
  }else
    refid <- -1

  if (!is.null(select)) {
    res <- filter(res, group_id %in% select)#select the target groups
  }else
  {
    if(refid>0)
      res <- filter(res, group_id != refid)#exclude the reference group if ref id is supplied

  }
  res <- group_by(res, group_id)

  res <- res %>%
    group_split() %>%#process each group
    map(function(df) {
      data <- df[[type]]#grab the key column
      if(length(data) == 1 && data == "")
        data <- character()
      unknown <- setdiff(data, ref) #find the item that are in the target groups but not in ref
      missing <- setdiff(ref, data) # find the item that are in the ref but not in the other group
      if (length(unknown) > 0 || length(missing) > 0) {
        list(unknown = unknown, missing = missing)
      }
    }) %>%
    set_names(group_keys(res)[["group_id"]])

  class(res) <- c(paste0("cqc_match_result_", type), "cqc_match_result", class(res))

  attr(res, "groups") <- x
  #perform the auto matching between the unknown and missing entries
  solution <- cqc_recommend(res, ...)

  res <- list(solution = solution, match_result = res, ref = ref)
  class(res) <- c("cqc_match_result_and_solution", class(res))
  res
}

#' Manual update/delete the match report
#'
#' There are cases that the automatic string match from `cqc_match` can't resolve.
#' This function provides the remedy to manually update the match result.
#'
#' cqc_match_update() matches the items from target groups to the reference by the named vector through 'map' argument.
#' cqc_match_remove() undo the matching that are produced either from cqc_match_update or cqc_match()
#' cqc_match_delete_unmatched() forces the unmatched items marked for deletion.
#'
#' @param x 'cqc_match_result_and_solution' object generated by "cqc_match"
#' @param group_id the group to update. If NULL, then operate on all groups
#' @param map For cqc_match_update, it is a named vector that provides the paired strings(i.e. old value vs new value).
#'            For cqc_match_remove or cqc_match_delete_unmatched, it is a character vector that provides the items to be removed from matched pairs or mark
#'            the unmatched items for deletion.
#' @examples
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' qc_cf_list <- cqc_load_fcs(fcs_files)
#' groups <- cqc_check(qc_cf_list, type = "marker")
#' match_result <- cqc_match(groups, ref = 3)
#'
#' # add the matching pairs to the auto-matching results
#' match_result <- cqc_match_update(match_result, map = c("PTPRC PE" = "CD45 PE"), group_id = 5)
#'
#' # You can delete any match added manually or automatically by cqc_match
#' match_result <- cqc_match_remove(match_result, "PTPRC PE", group_id = 5)
#'
#' # or mark selected items for deletion
#' #match_result <- cqc_drop_unmatched(match_result, "PTPRC PE", group_id = 5)
#'
#' # You can also omit the group_id to apply to all groups
#' match_result <- cqc_match_update(match_result, map = c("PTPRC PE" = "CD45 PE"))
#' @export
#' @rdname cqc_match_update
cqc_match_update <- function(x, map, group_id = NULL){
  if(!is(x, "cqc_match_result_and_solution"))
    stop("x must be the match result generated by cqc_match call!")
  if(length(names(map)) == 0)
    stop("'map' should be named character!")
  match_res <- x$match_result
  grp.attr <- attr(x$solution, "group")
  cls <- class(x$solution)
  type <- sub("cqc_solution_", "", cls[1])
  #if groupid is not given then apply to all groups
  if(is.null(group_id))
  {
    group_id <- names(match_res)
  }else{
    group_id <- as.character(group_id)
  }
  from <- names(map)#grab the src from new input

  #loop through each group
  new_solution <- map_dfr(group_id, function(gid){
    this_match_res <- match_res[[gid]]
    keys <- filter(grp.attr, group_id == gid)[[type]]
    thisgrp <- filter(x$solution, group_id == gid)
    existing.from <- thisgrp[["from"]]
    existing.to <- filter(thisgrp, !is.na(to))[["to"]]

    #validity check from
    idx <- from %in% x$ref
    if(any(idx))
      stop(paste(from[idx], collapse = ", "), " are references!Something is wrong with the 'map'! ")

    idx <- from %in% existing.from
    if(any(idx))
    {
      print(filter(thisgrp, from %in% !!from))
      stop("Found the existing match. Please call 'cqc_match_remove' to clear the existing match first")
    }
    #only add the entries that are in unknown list
    idx <- from %in% this_match_res[["unknown"]]
    from <- from[idx]
    to <- map[from]

    #validity check for to
    #skip to = NA, which is for deletion
    to.non_na <- to[!is.na(to)]
    if(length(is.na(to)) > 0)
    {
      idx <- to.non_na %in% x$ref
      if(any(!idx))
        stop(paste(to.non_na[!idx], collapse = ", "), " are not valid references!")
      idx <- to.non_na %in% existing.to
      if(any(idx))
      {
        print(filter(thisgrp, to %in% !!to))
        stop("Found the existing match. Please call 'cqc_match_remove' to clear the existing match first")
      }
      idx <- to.non_na %in% match_res[[gid]][["missing"]]
      if(any(!idx))
        stop(paste(to.non_na[!idx], collapse = ", "), " are already perfectly matched!")

    }
    #add the new pair into solution
    tb <- tibble(group_id = as.integer(gid), from = from, to = as.character(to))

    type <- sub("cqc_solution_", "", class(x$solution))[1]

    ## add missing items to insertion list (only for keywords task at the moment)
    if(type == "keyword")
    {
      #check if unknown list will be cleared with the new solution
      from.all <- c(existing.from, from)
      from.all <- from.all[!is.na(from.all)]
      if (setequal(from.all, this_match_res[["unknown"]])) {
          #mark the rest of un-matched missing ref as insertion
        todel <- this_match_res[["missing"]]
        #exclude the matched solution item
        to.all <- c(existing.to, to)
        to.all <- to.all[!is.na(to.all)]
        todel <- todel[!todel %in% to.all]
        if(length(todel) > 0)
          tb <- add_row(tb, group_id = as.integer(gid), from = NA, to = todel)

      }
    }
    tb
  })
  if(nrow(new_solution)==0)
    stop(paste(from, collapse = ", "), " not found!")
  #append the new solution to the existing one
  x$solution <- bind_rows(x$solution, new_solution)
  #restore attr
  attr(x$solution, "group") <- grp.attr
  class(x$solution) <- cls
  x
}

#' @export
#' @rdname cqc_match_update
cqc_match_delete_unmatched <- function(x, map, group_id = NULL){
  todel <- rep(NA, length(map))
  names(todel) <- map
  cqc_match_update(x, todel)

}
#' @export
#' @rdname cqc_match_update
cqc_match_remove <- function(x, map, group_id = NULL){
  if(!is(x, "cqc_match_result_and_solution"))
    stop("x must be the match result generated by cqc_match call!")
  if(length(names(map)) > 0)
    stop("'map' should be unnamed character!")
  match_res <- x$match_result
  grp.attr <- attr(x$solution, "group")
  cls <- class(x$solution)

  exclude <- filter(x$solution, from %in% map)
  if(!is.null(group_id))
  {
    exclude <- filter(exclude, group_id %in% group_id)
  }
  if(nrow(exclude) > 0)
  {
    x$solution <- anti_join(x$solution, exclude, by = c("group_id", "from"))
    #restore attr
    attr(x$solution, "group") <- grp.attr
    class(x$solution) <- cls
  }else
    stop("No existing matches for ", paste(map, collapse = ", "), " to delete!")
  x
}
#' #' Convert the match_result to table
#' #' @param x "cqc_match_result" returned as part of cqc_match call
#' #' @importFrom dplyr as_tibble
#' #' @export
#' as_tibble.cqc_match_result <- function(x, ...) {
#'   map_dfr(x, function(i) {
#'     tibble(
#'       "Not in reference" = paste(i[["unknown"]], collapse = ","),
#'       "Missing channels" = paste(i[["missing"]], collapse = ",")
#'     )
#'   }, .id = "group_id")
#' }
