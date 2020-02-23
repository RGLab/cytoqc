context("cqc_check")
library(flowCore)
data(GvHD)
cflist <- fsApply(GvHD[1:2], function(fr)flowFrame_to_cytoframe(fr[1:2, c(1:3,6)]), simplify = FALSE)
cqc_data <- cqc_cf_list(cflist)

test_that("cqc_fix_panel -- rotated", {
  #fix rotated chnl
  cf <- cqc_data[[1]]
  su <- range(cf, "data")
  cn <- colnames(cf)

  newcn <- cn[c(2,1,3,4)]
  newcn[3] <- "D"
  colnames(cf) <- newcn
  cf <- cf[,4:1]

  groups <- cqc_check(cqc_data, "panel")
  expect_equal(ncol(format(groups, "marker")), 3)

  cqc_fix_panel(groups, 2, "marker")
  expect_equal(colnames(cf),  cn[4:1])
  expect_equal(colnames(cqc_data[[2]]),  cn)
  #make sure the underlying data is correct as well
  expect_equal(range(cf, "data")[,colnames(su)], su)
  #fix rotated marker
  # pnl <- cf_get_panel(cqc_data[[1]])
  # cn <- pnl[["marker"]]
  # newcn <- cn[c(2,3,1)]
  # names(newcn) <- pnl[["channel"]]
  # markernames(cqc_data[[1]]) <- newcn
  # groups <- cqc_check(cqc_data, "panel")
  # expect_equal(ncol(format(groups, "marker")), 3)
  #
  # cqc_fix_panel(groups, 1, "marker")
  # expect_equal(colnames(cqc_data[[1]]),  newcn)

})
