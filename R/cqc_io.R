#' Load FCS files
#' load fcs into 'cqc_cf_list' object which is a list of cytoframes
#' @import flowWorkspace
#' @export
cqc_load_fcs <- function(files, is_h5 = TRUE, ...){
  res <- sapply(files, function(file)load_cytoframe_from_fcs(file, is_h5 = is_h5, ...))
  names(res) <- basename(names(res))
  attr(res, "class") <- "cqc_cf_list"
  res
}
#' @export
write_fcs <- function(x, ...)UseMethod("write_fcs")

#' The helper function to write the cleaned cqc_cf_list back to fcs
#' @param x cqc_cf_list
#' @param out the output directory that the FCS will be written
#' @export
write_fcs.cqc_cf_list <- function(x, out, verbose = TRUE){
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
