% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cqc_fix.R
\name{cqc_fix_panel}
\alias{cqc_fix_panel}
\title{Update panel info}
\usage{
cqc_fix_panel(x, ref, by, select = NULL)
}
\arguments{
\item{x}{"cqc_check_panel" result from 'cqc_check_panel' or 'cqc_check'call}

\item{ref}{the reference group id or a custom panel that is represented as tibble contains channel and marker columns (typically returned by 'cf_get_panel` )}

\item{by}{the column used as the fixed anchor, either 'channel' or 'marker'}

\item{select}{the groups to be updated. Default is all}
}
\description{
It updates either 'channel' or 'marker' info depends on the value of 'ref.col'
if ref.col is set to 'channel', then 'channel' column from 'panel' is used as the
reference to be matched against 'cqc_data' and 'marker'  will be set according to 'panel'
}
\examples{
\dontrun{
     res <- cqc_check(cqc_data, "panel")
     cqc_fix_panel(res, 1, "channel") #get panel info from group 1 and set markers of other groups  by using channel as reference

}
}
