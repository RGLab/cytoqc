
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

      #(1st mat)
      # Levenshtein (edit) distance
      dist_mat <- adist(refs, targets, ignore.case = TRUE)

      #(2nd mat)
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
      # (3rd mat)
      #check if one if the substring of the other
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

      #combine the info from 3 matrices above to find the best match

      #first scan the substr mat to get pairs that has substring match
      is_sub_idx <- which(is_substring_mat)
      #order these substring matched pairs by the approximate string dist
      #to break tie when multiples are substr-matched to one
      is_sub_idx <- is_sub_idx[order(dist_mat[is_sub_idx])]
      #then process the pairs that has no substr relations
      #and also discard the pairs that do not pass dist threshold
      no_sub_idx <- which((!is_substring_mat)&is_pass_mat)
      #order them by dist
      no_sub_idx <- no_sub_idx[order(dist_mat[no_sub_idx])]
      #add the pairs selected/ordered by the above rules
      for(idx in c(is_sub_idx, no_sub_idx))
      {
        if (length(refs_queue) == 0 || length(targets_queue) == 0)
          break#terminate the loop if one of the queues is empty

        # try to parse out the names of the pair by its x, y coordinates within the mat
        ridx <- idx %% nrows
        if (ridx == 0) {
          ridx <- nrows
        }
        cidx <- ceiling(idx / nrows)
        # get the pair
        from <- targets[cidx]
        to <- refs[ridx]
        #check if they are already matched by previous iterations
        tind <- match(from, targets_queue)
        rind <- match(to, refs_queue)
        #if new match then add them to the solution
        if(!is.na(tind) && !is.na(rind) )
        {
          df <- add_row(df, from = from, to = to)

          # pop the processed item
          refs_queue <- refs_queue[-rind]
          targets_queue <- targets_queue[-tind]
        }

      }#end for
    }
    #add those extra items to deletion list
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

#' @importFrom tibble as_tibble
#' @importFrom dplyr select distinct 
#' @importFrom tidyr pivot_wider replace_na
cluster_panels <- function(res, missing_penalty = 1.0){
  clipped <- as_tibble(res) %>%
    select(-object) %>%
    distinct() %>%
    # arrange(group_id, channel)
    pivot_wider(names_from = channel, values_from=marker)
    
  df <- t(as.data.frame(clipped %>% select(-c(group_id, nObject))))
  colnames(df) <- paste0("Group ", clipped$group_id)
  # Determine max distance for normalization
  max_dist <- max(replace_na(adist(df), 0))
  
  # Distance metric -- partial allows for substring match, which is asymmetric
  dm <- function(x,y){sum(replace_na(diag(adist(x, y, partial=TRUE)), missing_penalty*max_dist))}
  
  ngroups <- ncol(df)
  dist_mat <- matrix(apply(expand.grid(1:ngroups,1:ngroups), 1, function(pair){dm(df[,pair[[1]]], df[,pair[[2]]])}), nrow = ngroups)
  
  # To handle the asymmetry and give credit for substring match while still penalizing mismatch,
  # average across diagonal -- could add weight factors to give more credit to substring match
  upper <- lower <- dist_mat
  lower[upper.tri(dist_mat)] <- 0
  upper[lower.tri(dist_mat)] <- 0 
  dist_mat <- (t(upper) + lower) / 2
  
  # normalize by largest value to make height simple [0,1]
  max_dist <- max(dist_mat)
  if(max_dist > 0)
    dist_mat <- dist_mat / max_dist
  
  colnames(dist_mat) <- rownames(dist_mat) <- colnames(df)
  dist_mat <- as.dist(dist_mat)
  clustered <- hclust(dist_mat)
}
