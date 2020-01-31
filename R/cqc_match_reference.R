#' find the the difference between the reference and target group
#'
#'
#' @param x cqc_prepare object returned by cqc_prepare call
#' @param ...
#'        ref specifies the reference, which can be either an integer group id or a characte vector giving the actual values of the reference
#'        select the group ids selected for processing
#'        type the qc type (either "channle", "marker", "gate")
#' @export
cqc_match <- function(x, ...) UseMethod("cqc_match")

#' @export
cqc_match.cqc_prepare_channel <- function(x, ...) {
  res <- match_reference(x, type = "channel", ...)
  class(res) <- c("cqc_match_result_channel", class(res))
  res
}

#' @export
cqc_match.cqc_prepare_marker <- function(x, ...) {
  res <- match_reference(x, type = "marker", ...)
  class(res) <- c("cqc_match_result_marker", class(res))
  res
}

#' @export
cqc_match.cqc_prepare_keyword <- function(x, ...) {
  res <- match_reference(x, type = "keyword", ...)
  class(res) <- c("cqc_match_result_keyword", class(res))
  res
}
#' @export
cqc_match.cqc_prepare_gate <- function(x, ...) {
  res <- match_reference(x, type = "gate", ...)
  class(res) <- c("cqc_match_result_gate", class(res))
  res
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
match_reference <- function(x, ref, select = NULL, type, delimiter = "|") {
  res <- summary(x)
  if (is.numeric(ref)) {
    ref <- res %>%
      filter(group_id %in% ref) %>%
      pull(type)
  }

  if (!is.null(select)) {
    res <- filter(res, group_id %in% select)
  }
  res <- group_by(res, group_id)

  res <- res %>%
    group_split() %>%
    map(function(df) {
      data <- df[[type]]
      unknown <- setdiff(data, ref)
      missing <- setdiff(ref, data)
      if (length(unknown) > 0 || length(missing) > 0) {
        list(unknown = unknown, missing = missing)
      }
    }) %>%
    set_names(group_keys(res)[["group_id"]])

  class(res) <- c("cqc_match_result", class(res))
  attr(res, "groups") <- x

  res
}


#' @importFrom dplyr as_tibble
#' @export
as_tibble.cqc_match_result <- function(x) {
  map_dfr(x, function(i) {
    tibble(
      "Not in reference" = paste(i[["unknown"]], collapse = ","),
      "Missing channels" = paste(i[["missing"]], collapse = ",")
    )
  }, .id = "group_id")
}
