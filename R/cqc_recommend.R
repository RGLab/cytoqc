
#' @export
cqc_recommend <- function(x, ...) UseMethod("cqc_recommend")
#' @export
cqc_recommend.cqc_check <- function(x, max.distance = 0.1, ...) {
  match_result <- cqc_match(groups, ...)
  cqc_recommend(match_result, max.distance = max.distance)
}

#' @export
cqc_recommend.cqc_match_result_channel <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_channel", attr(res, "class"))
  res
}
#' @export
cqc_recommend.cqc_match_result_marker <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_marker", attr(res, "class"))
  res
}
#' @export
cqc_recommend.cqc_match_result_keyword <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_keyword", attr(res, "class"))
  res
}
#' @export
cqc_recommend.cqc_match_result_gate <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_gate", attr(res, "class"))
  res
}
#' Find solution to resolve the discrepancy discovered by match_reference
#'
#' It tries to find the aproximate match(based on 'agrep') between the target and reference as well as the extra redundunt items that can be removed.
#' @name cqc_recommend
#' @return a table (with 'from' and 'to' columns) represents the itemized fix recommendation. When 'to' is 'NA', it means the entry is redundunt and can be removed
#' @examples
#' \dontrun{
#' solution <- cqc_recommend(groups, select = c(1, 4))
#' }
#' @param x A CQC object of some kind. See vignettes.
#' @param max.distance Maximum distance allowed for a match. See ?agrep
#' @param ... additional arguments not for the user.
#' @importFrom tibble tibble add_row
#' @importFrom utils adist
#' @export
cqc_recommend.cqc_match_result <- function(x, max.distance = 0.1, ...) {
  res <- map_dfr(x, function(check_result) {
    unknown <- check_result[["unknown"]]
    missing <- check_result[["missing"]]
    df <- tibble(from = character(), to = character())
    # iteratively try to find the aproximate match between two vecs
    while (length(unknown) > 0) {
      if (length(missing) > 0) {
        # Levenshtein (edit) distance
        dist_mat <- adist(unknown, missing, ignore.case = TRUE)
        nrows <- nrow(dist_mat)
        ncols <- ncol(dist_mat)
        # pick the best match
        idx <- which.min(dist_mat)
        # browser()
        # get x, y coordinates
        ridx <- idx %% nrows
        if (ridx == 0) {
          ridx <- nrows
        }
        cidx <- ceiling(idx / nrows)
        # get the pair
        from <- unknown[ridx]
        to <- missing[cidx]
        # pop the processed item
        unknown <- unknown[-ridx]
        missing <- missing[-cidx]
        # check if exceeds max.distance
        # agrep can be avoided if the formula of max.distance used by agrep is figured out
        ind <- agrep(from, to, ignore.case = TRUE, max.distance = max.distance)
        if (length(ind) > 0) { #
          df <- add_row(df, from = from, to = to)
        } else {
          next
        } # otherwise stop the matching process
      } else # extra params
      {
        df <- add_row(df, from = unknown, to = NA)
        unknown <- character()
      }
    }
    df
  }, .id = "group_id") %>% mutate(group_id = as.integer(group_id))
  attr(res, "class") <- c("cqc_solution", attr(res, "class"))
  attr(res, "group") <- attr(x, "group")
  res
}
