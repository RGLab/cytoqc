context("cqc_check")
library(flowCore)
data(GvHD)

test_that("cqc_fix_panel -- rotated", {
  cflist <- fsApply(GvHD[1:2], function(fr)flowFrame_to_cytoframe(fr[1:2, c(1:3,6)]), simplify = FALSE)
  cqc_data <- cqc_cf_list(cflist)
  #fix rotated chnl
  cf <- cqc_data[[1]]
  su <- range(cf, "data")
  cn <- colnames(cf)

  newcn <- cn[c(2,1,3,4)]
  newcn[3] <- "D"
  colnames(cf) <- newcn
  cf <- cf[,4:1]
  #craft a spillover matrix
  sp <- matrix(rnorm(9), nrow = 3, dimnames = list(NULL, newcn[1:3]))
  keyword(cf)[["SPILL"]] <- sp
  groups <- cqc_check(cqc_data, "panel")
  expect_equal(ncol(format(groups, "marker")), 3)

  cqc_fix_panel(groups, 2, "marker")
  expect_equal(colnames(cf),  cn[4:1])
  expect_equal(colnames(cqc_data[[2]]),  cn)
  #make sure the underlying data is correct as well
  expect_equal(range(cf, "data")[,colnames(su)], su)
  expect_equal(colnames(keyword(cf)[["SPILL"]]), cn[1:3])

  #fix rotated marker
  pnl <- cf_get_panel(cqc_data[[1]])
  marker <- pnl[["marker"]]
  newmarker <- marker[c(2,3,1)]
  names(newmarker) <- pnl[["channel"]][1:3]
  markernames(cqc_data[[1]]) <- newmarker
  groups <- cqc_check(cqc_data, "panel")
  expect_equal(ncol(format(groups)), 3)

  cqc_fix_panel(groups, 1, "marker")
  pnl <- cf_get_panel(cqc_data[[1]])
  expect_equal(pnl,  cf_get_panel(cqc_data[[2]]))


  })
