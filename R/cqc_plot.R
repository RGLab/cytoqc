#' Specialized plots for objects from the \code{cytoqc} package
#' 
#' This is a generic function for plotting functionality specific
#' to objects from the \code{\link[=cytoqc-package]{cytoqc}} package.
#' 
#' It is currently only implemented for results of \code{\link{cqc_cluster}}.
#' 
#' @param y object to be plotted
#' @param ... additional arguments (currently unused)
#' @export
cqc_plot <- function(y, ...) UseMethod("cqc_plot")

#' @param y an object of class cqc_cluster resulting from a cqc_cluster call
#' @rdname cqc_cluster
#' @export
cqc_plot.cqc_cluster <- function(y, ...){
  # Maybe look in to how to override stats:::plot.hclust
  # or just work with the args there
  plot(y$cluster_results, 
       xlab="Group Number",
       sub="")
  if(!is.null(y$height))
    abline(a=y$height, b=0, col="blue", lty=2)
  if(y$ngroup > 1 && y$ngroup < length(unique(y$group_membership$old_group)))
    rect.hclust(y$cluster_results, k = y$ngroup, cluster=y$plot_groups, border = "red")
}