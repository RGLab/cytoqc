#' Load cqc_data
#' 
#' load fcs or h5 files into \code{\link{cqc_cf_list}} object which is a list of \code{\link[flowWorkspace]{cytoframe}} objects.
#' This is the method to construct the core data object for \code{\link[cytoqc:cytoqc-package]{cytoqc}}.
#' @param files the fcs or h5 file paths
#' @param is_h5 \code{logical} should the cytoframe be constructed as an h5 disk-backed structure. Default \code{TRUE}.
#'                              It is ignored for \code{cqc_load_h5}
#' @param ... parameters passed to 'load_cytoframe_from_fcs' or 'load_cytoframe_from_h5'
#' @rdname cqc_load_fcs
#' @import flowWorkspace
#' @importFrom methods is
#' @examples 
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' cqc_cf_list <- cqc_load_fcs(fcs_files)
#' @export
cqc_load_fcs <- function(files, is_h5 = TRUE, ...) {
  res <- sapply(files, function(file) load_cytoframe_from_fcs(file, is_h5 = is_h5, ...))
  names(res) <- basename(names(res))
  cqc_cf_list(res)
}

#' @rdname cqc_load_fcs
#' @examples
#' \dontrun{
#' h5_files <- list.files(system.file("extdata", "gs_bcell_auto", package = "flowWorkspaceData"),
#'                        pattern = ".h5", full.names = TRUE)
#' cqc_cf_list <- cqc_load_h5(h5_files)
#' }
#' @export
cqc_load_h5 <- function(files, is_h5 = TRUE, ...) {
  res <- sapply(files, function(file) load_cytoframe_from_h5(file, readonly = FALSE, ...))
  names(res) <- basename(names(res))
  cqc_cf_list(res)
}

#' Construct a \code{cqc_cf_list} object from a list of \code{cytoframe} objects
#'
#' This is the core data object for \code{\link[cytoqc:cytoqc-package]{cytoqc}}.
#'
#' @param x a named list of \code{\link[flowWorkspace]{cytoframe}} objects
#' @examples 
#' # This is just for illustration. cqc_load_fcs will normally take care of this step.
#' fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
#' cf_list <- lapply(fcs_files[1:3], load_cytoframe_from_fcs)
#' names(cf_list) <- fcs_files[1:3]
#' 
#' # Construct a cqc_cf_list object from a list of cytoframes
#' cf_list <- cqc_cf_list(cf_list)
#' 
#' @export
cqc_cf_list <- function(x) {
  if (!is.list(x)) {
    stop("x must be a list!")
  }
  if (length(names(x)) != length(x)) {
    stop("x must be a named list!")
  }

  for (i in x) {
    if (!is(i, "cytoframe")) {
      stop("Each element in x must be a cytoframe!")
    }
  }
  class(x) <- c("cqc_cf_list", "cqc_data", class(x))
  x
}

#' Construct a 'cqc_gs' object from a 'GatingSet'
#'
#' This is mainly for distpatch and wrapping many of the
#' qc operations on a \code{\link{cqc_cf_list}} object for the GatingSet's
#' underlying data
#' @param x a GatingSet object
#' @examples
#' gs_path <- system.file("extdata", "gslist_manual_QC", "gs1", package = "cytoqc")
#' gs <- load_gs(gs_path)
#' qc_gs <- cqc_gs(gs)
#' @export
cqc_gs <- function(x) {
  if (!is(x, "GatingSet")){
    stop("x must be a GatingSet!")
  }
  # Expand GatingSet in to a list of GatingHierarchies for gate qc
  # Have to avoid a get_cytoset call which could error with internal inconsistency in channel
  ghlist <- lapply(1:length(flowWorkspace:::.cpp_getSamples(x@pointer)), function(idx) {x[[idx]]})
  names(ghlist) <- sampleNames(x)

  class(ghlist) <- c("cqc_gs", "cqc_data", class(ghlist))
  ghlist
}


#' Write out tidied flow data (\code{cqc_cf_list}) back to fcs or h5
#' @param x cqc_cf_list
#' @param out the output directory that the FCS or h5 will be written
#' @param verbose whether to print each sample name during the writing process
#' @param ... other arguments passed down to 'write.FCS'
#' @importFrom flowCore write.FCS
#' @rdname cqc_write_fcs
#' @export
cqc_write_fcs <- function(x, out, verbose = TRUE, ...) {
  if (!dir.exists(out)) {
    dir.create(out)
  }
  for (sn in names(x))
  {
    if (verbose) {
      message("writing ", sn)
    }
    fr <- x[[sn]]
    write.FCS(fr, filename = file.path(out, sn), ...)
  }
}

#' @rdname cqc_write_fcs
#' @export
cqc_write_h5 <- function(x, out, verbose = TRUE) {
  if (!dir.exists(out)) {
    dir.create(out)
  }
  for (sn in names(x))
  {
    if (verbose) {
      message("writing ", sn)
    }
    fr <- x[[sn]]
    cf_write_h5(fr, filename = file.path(out, sn))
  }
}


#' Construct a 'cqc_gs_list' object from a list of 'GatingSet' objects
#'
#' For the methods dispatching purpose
#'
#' @param x a list of 'GatingSet' objects
#' @examples 
#' gs_paths <- list.files(system.file("extdata", "gslist_manual_QC", package = "cytoqc"), full.names = TRUE)
#' gs1 <- load_gs(gs_paths[[1]])
#' gs2 <- load_gs(gs_paths[[2]])
#' qc_gs_list <- cqc_gs_list(list(gs1, gs2))
#' groups <- cqc_check(qc_gslist, type="gate")
#' 
#' @export
cqc_gs_list <- function(x) {
  if (!is.list(x)) {
    stop("x must be a list!")
  }
  # browser()
  names(x) <- lapply(x, function(gs) {
    if (!is(gs, "GatingSet")) {
      stop("Each element in x must be a GatingSet!")
    }
    identifier(gs)
  })
  class(x) <- c("cqc_gs_list", "cqc_data", class(x))
  x
}
