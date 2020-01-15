#' Load FCS files
#' load fcs into 'cqc_cf_list' object which is a list of cytoframes.
#' This is the method to construct the core data object for cytoQC.
#' @param is_h5 \code{logical} should the cytoframe be constructed as an h5 disk-backed structure. Default \code{TRUE}.
#' @import flowWorkspace
#' @importFrom methods is
#' @export
cqc_load_fcs <- function(files, is_h5 = TRUE, ...){
  res <- sapply(files, function(file)load_cytoframe_from_fcs(file, is_h5 = is_h5, ...))
  names(res) <- basename(names(res))
  cqc_cf_list(res)
}

#' Construct a 'cqc_cf_list' object from a list of 'cytoframe' objects
#'
#' This is the core data object for CytoQC.
#'
#' @param x a list of 'cytoframe' objects
#' @export
cqc_cf_list <- function(x){

  if(!is.list(x))
    stop("x must be a list!")
  if(length(names(x))!=length(x))
    stop("x must be a named list!")

  for(i in x)
    if(!is(i, "cytoframe"))
      stop("Each element in x must be a cytoframe!")
  class(x) <- c("cqc_cf_list", class(x))
  x

}

#' @export
write_fcs <- function(x, ...)UseMethod("write_fcs")

#' Write out tidied flow data (cqc_cf_list) back to fcs
#' @param x cqc_cf_list
#' @param ... additional arguments.
#' @param out the output directory that the FCS will be written
#' @importFrom flowCore write.FCS
#' @export
write_fcs.cqc_cf_list <- function(x, out, verbose = TRUE,...){
  if(!dir.exists(out))
    dir.create(out)
  for(sn in names(x))
  {
    if(verbose)
      message("writing ", sn)
    # fr <- cytoframe_to_flowFrame(x[[sn]])
    fr <- x[[sn]]
    write.FCS(fr, filename = file.path(out, sn))
  }

}

#' Construct a 'cqc_gs_list' object from a list of 'GatingSet' objects
#'
#' For the methods dispatching purpose
#'
#' @param x a list of 'GatingSet' objects
#' @export
cqc_gs_list <- function(x){

  if(!is.list(x))
    stop("x must be a list!")
  # browser()
  names(x) <- lapply(x, function(gs){
                if(!is(gs, "GatingSet"))
                  stop("Each element in x must be a GatingSet!")
                identifier(gs)
              })
  class(x) <- c("cqc_gs_list", class(x))
  x

}
