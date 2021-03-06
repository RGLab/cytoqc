% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cqc_check.R
\name{cqc_get_data}
\alias{cqc_get_data}
\title{Extract the data from the result of a \code{cqc_check} call.}
\usage{
cqc_get_data(groups, id = NULL)
}
\arguments{
\item{groups}{the object returned by \code{\link{cqc_checks}}}

\item{id}{the group id to be selected from the dataset, default is NULL, meaning all data}
}
\description{
Extract the data from the result of a \code{cqc_check} call.
}
\examples{
fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
qc_cf_list <- cqc_load_fcs(fcs_files)
channel_groups <- cqc_check(qc_cf_list, type = "channel")
summary(channel_groups)
group_3_cf_list <- cqc_get_data(channel_groups, 3)
# A list of cytoframes
group_3_cf_list
}
