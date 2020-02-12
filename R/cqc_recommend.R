
cqc_recommend <- function(x, ...) UseMethod("cqc_recommend")
cqc_recommend.cqc_check <- function(x, max.distance = 0.1, ...) {
  match_result <- cqc_match(groups, ...)
  cqc_recommend(match_result, max.distance = max.distance)
}

cqc_recommend.cqc_match_result_channel <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_channel", attr(res, "class"))
  res
}
cqc_recommend.cqc_match_result_marker <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_marker", attr(res, "class"))
  res
}
cqc_recommend.cqc_match_result_keyword <- function(x, ...) {
  res <- cqc_recommend.cqc_match_result(x, ...)
  attr(res, "class") <- c("cqc_solution_keyword", attr(res, "class"))
  res
}
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
#' @param partial whether to do the partial sub string matching before the approximate string matching
#' @param ... additional arguments not for the user.
#' @importFrom tibble tibble add_row
#' @importFrom utils adist
cqc_recommend.cqc_match_result <- function(x, max.distance = 0.1, partial = TRUE, ...) {
  res <- map_dfr(x, function(check_result) {
    targets_queue <- targets <- check_result[["unknown"]]
    refs_queue <- refs <- check_result[["missing"]]
    df <- tibble(from = character(), to = character())

    if (length(refs_queue) > 0 && length(targets_queue) > 0)
    {


      # Levenshtein (edit) distance
      dist_mat <- adist(refs, targets, ignore.case = TRUE)

      # check if each match passes the max.distance check
      # agrep can be avoided if the formula of max.distance used by agrep is figured out
      is_pass_mat <- do.call(rbind, sapply(refs, function(ref){
                                     sapply(targets, function(target){
                                        #default definition for max.distance of agrep is non-symetric
                                        is.match <- agrepl(ref, target, ignore.case = TRUE, max.distance = max.distance)
                                        is.match || agrepl(target, ref, ignore.case = TRUE, max.distance = max.distance)

                                       })
                                    }, simplify = FALSE)
                            )

      #check if substring
      if(partial)
      {
        is_substring_mat <- do.call(rbind, sapply(refs, function(ref){
          sapply(targets, function(target){
            grepl(paste0("\\Q", ref, "\\E"), target, ignore.case = TRUE)||grepl(paste0("\\Q", target, "\\E"), ref, ignore.case = TRUE)

          })
        }, simplify = FALSE)
        )
      }else#create dummy mat
        is_substring_mat <- matrix(TRUE, nrow = length(refs), ncol = length(targets), dimnames = list(refs, targets))

      nrows <- nrow(dist_mat)
      ncols <- ncol(dist_mat)


      #first scan the substr matched entries
      is_sub_idx <- which(is_substring_mat)
      #order by dist
      is_sub_idx <- is_sub_idx[order(dist_mat[is_sub_idx])]
      #then scan the non-substr  entries
      no_sub_idx <- which((!is_substring_mat)&is_pass_mat)#filter by dist threshold mat
      no_sub_idx <- no_sub_idx[order(dist_mat[no_sub_idx])]
      for(idx in c(is_sub_idx, no_sub_idx))
      {
        if (length(refs_queue) == 0 || length(targets_queue) == 0)
          break#terminate the loop if one of the queues is empty

          # get x, y coordinates
        ridx <- idx %% nrows
        if (ridx == 0) {
          ridx <- nrows
        }
        cidx <- ceiling(idx / nrows)
        # get the pair
        from <- targets[cidx]
        to <- refs[ridx]
        tind <- match(from, targets_queue)
        rind <- match(to, refs_queue)

        if(!is.na(tind) && !is.na(rind) )#skip already-matched pairs
        {
          df <- add_row(df, from = from, to = to)

          # pop the processed item
          refs_queue <- refs_queue[-rind]
          targets_queue <- targets_queue[-tind]
        }

      }#end for
    }

    # iteratively try to find the aproximate match between two vecs
    if (length(targets_queue) > 0 && length(refs_queue) == 0) {

         df <- add_row(df, from = targets_queue, to = NA)
         targets_queue <- character()

    }
    df
  }, .id = "group_id") %>% mutate(group_id = as.integer(group_id))
  attr(res, "class") <- c("cqc_solution", attr(res, "class"))
  attr(res, "group") <- attr(x, "group")
  res
}
