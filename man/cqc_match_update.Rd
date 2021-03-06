% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cqc_match.R
\name{cqc_match_update}
\alias{cqc_match_update}
\alias{cqc_match_delete_unmatched}
\alias{cqc_match_remove}
\title{Manual update/delete the match report}
\usage{
cqc_match_update(x, map, group_id = NULL)

cqc_match_delete_unmatched(x, map, group_id = NULL)

cqc_match_remove(x, map, group_id = NULL)
}
\arguments{
\item{x}{'cqc_match_result_and_solution' object generated by "cqc_match"}

\item{map}{For cqc_match_update, it is a named vector that provides the paired strings(i.e. old value vs new value).
For cqc_match_remove or cqc_match_delete_unmatched, it is a character vector that provides the items to be removed from matched pairs or mark
the unmatched items for deletion.}

\item{group_id}{the group to update. If NULL, then operate on all groups}
}
\description{
There are cases that the automatic string match from `cqc_match` can't resolve.
This function provides the remedy to manually update the match result.
}
\details{
cqc_match_update() matches the items from target groups to the reference by the named vector through 'map' argument.
cqc_match_remove() undo the matching that are produced either from cqc_match_update or cqc_match()
cqc_match_delete_unmatched() forces the unmatched items marked for deletion.
}
\examples{
fcs_files <- list.files(system.file("extdata", "GvHD_QC", package = "cytoqc"), full.names = TRUE)
qc_cf_list <- cqc_load_fcs(fcs_files)
groups <- cqc_check(qc_cf_list, type = "marker")
match_result <- cqc_match(groups, ref = 3)

# add the matching pairs to the auto-matching results
match_result <- cqc_match_update(match_result, map = c("PTPRC PE" = "CD45 PE"), group_id = 5)

# You can delete any match added manually or automatically by cqc_match
match_result <- cqc_match_remove(match_result, "PTPRC PE", group_id = 5)

# or mark selected items for deletion
#match_result <- cqc_drop_unmatched(match_result, "PTPRC PE", group_id = 5)

# You can also omit the group_id to apply to all groups
match_result <- cqc_match_update(match_result, map = c("PTPRC PE" = "CD45 PE"))
}
