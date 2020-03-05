#' Cluster groups based on the result of cqc_check
#' 
#' Sets of labels (strings representing channels, markers, keywords, or panels consisting of channel/marker pairs) for each group
#' will be clustered based on the summed distances of their aligned elements. This provides an automated way to look for groups with
#' sets of labels that likely represent the same underlying values but with minor differences or errors.
#'
#' @param x \code{cqc_check} result returned by cqc_check call
#' @param ...
#'        ngroup: integer specifying number of groups in to which the clusters will be combined
#'        
#'        height: number between 0 and 1 that specifies a cut height for the tree to determine the groups. Ignored if ngroup is specified.
#'        
#'        missing_penalty: multiplicative factor used to impose penalty for labels missing from a group. The distance for a missing label
#'        will be calculated as missing_penalty times the largest distance between non-missing labels. Defaults to 1.
#' @return a cqc_cluster object whose primary value is the tibble in its group_membership slot providing the map between the original
#' groups and the combined groups resulting from clustering
#' @export
cqc_cluster <- function(x, ...) UseMethod("cqc_cluster")

#' @export
cqc_cluster.default <- function(x, ...) {
  # Improve this error message
  stop("The input is not a valid 'cqc_check' result!\nPlease make sure to follow the right order of 'cqc_check-->cqc_cluster-->cqc_split' workflow!")
}

#' @export
cqc_cluster.cqc_check_channel <- function(x, ...){
  res <- cluster_labels(x, type = "channel", ...)
  res
}

#' @export
cqc_cluster.cqc_check_marker <- function(x, ...){
  res <- cluster_labels(x, type = "marker", ...)
  res
}

#' @export
cqc_cluster.cqc_check_keyword <- function(x, ...){
  res <- cluster_labels(x, type = "keyword", ...)
  res
}

#' @export
cqc_cluster.cqc_check_panel <- function(x, ...){
  res <- cluster_labels(x, type = "panel", ...)
  res
}



#' @importFrom tibble as_tibble
#' @importFrom dplyr select distinct 
#' @importFrom tidyr pivot_wider replace_na
#' @export
cluster_labels <- function(check_res, type, ngroup = NULL, height = NULL, missing_penalty = 1.0){
  if(type != "panel"){
    clipped <- summary(check_res) %>%
      group_by(group_id) %>%
      mutate(label_idx=1:n()) %>%
      ungroup() %>%
      pivot_wider(names_from = label_idx, values_from = !!type)
  }else{
    clipped <- check_res %>%
      select(-object) %>%
      distinct() %>%
      pivot_wider(names_from = channel, values_from=marker)
  }
  
  # hclust will succeed, but stats:::plot.hclust yields opaque error when called on result and it doesn't really make sense
  # to be here if you have 1 or 2 groups anyway
  if(nrow(clipped) < 3)
    stop(paste0("This cqc_check object only has ", nrow(clipped), " groups. Clustering of groups requires 3 or more groups."))
  
  df <- t(as.data.frame(clipped %>% select(-c(group_id, nObject))))
  colnames(df) <- clipped$group_id
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
  
  new_groups <- NULL
  if (!is.null(ngroup))
    new_groups <- cutree(clustered, k = ngroup)
  else if (!is.null(height))
    new_groups <- cutree(clustered, h = height)
  else
    new_groups <- cutree(clustered, k = ngroups)
  
  if(is.null(ngroup))
    ngroup <- length(unique(new_groups))
  
  group_membership <- tibble(old_group=as.integer(names(new_groups)), new_group=new_groups) %>%
    arrange(old_group)
  
  cluster_res <- list("cluster_results"=clustered, "group_membership"=group_membership, "ngroup"=ngroup, "height"=height, "plot_groups"=new_groups)
  class(cluster_res) <- c("cqc_cluster", class(cluster_res))
  class(cluster_res) <- c(paste0("cqc_cluster_", type), class(cluster_res))
  attr(cluster_res, "check_res") <- check_res
  cluster_res
}


#' Split the result of 'cqc_cluster' into groups
#'
#' It is used to split samples into separate groups when they can't be reconciled into the same group.
#'
#' @importFrom purrr walk
#' @param x cqc_cluster object
#' @param f,drop,... not used
#' @export
split.cqc_cluster <- function(x, f, drop = FALSE, ...) {
  check_res <- attr(x, "check_res")
  cqc_data <- attr(check_res, "data")
  data_type <- class(cqc_data)
  vec <- check_res %>%
    select(c(object, group_id)) %>%
    distinct() %>%
    rename(old_group=group_id) %>%
    inner_join(x$group_membership, by = "old_group") %>%
    pull(new_group)
  if (length(vec) != length(unique(check_res$object)))
    stop("Invalid cqc_cluster object. Not all samples have new group assignment.")
  split(cqc_data, vec) %>% map(function(i) {
    class(i) <- c(data_type, class(i))
    i
  })
}


