data(GvHD)
fs <- GvHD[1:2]
fr <- GvHD[[1]]

test_that("change colnames separately through parameter or exprs", {
  expect_error(pData(parameters(fr))[["name"]][[1]] <- "A", "don't match")
  expect_error(colnames(exprs(fr))[1] <- "A", "unable to replace")
  colnames(fr@exprs)[1] <- "A"
  expect_error(fr[, "FSC-H"], "out of bounds")
  expect_error(fr[, "A"], "out of bounds")
})

test_that("colnames", {
  channels <- c('FSC-H', 'SSC-H', 'FL1-H', 'FL2-H', 'FL3-H', 'FL2-A', 'FL4-H', 'Time')
  expect_equal(colnames(fr), channels)
  
  expect_equal(colnames(fs), channels)
  expect_equal(colnames(fr[, c(3,5)]), channels[c(3,5)])
  
})

test_that("colnames<-", {
  chnls <- c("A", "B")

  #update colnames for flowFrame
  colnames(fr)[c(1,3)] <- chnls
  expect_equal(colnames(fr)[c(1,3)], chnls)
  expect_equivalent(unlist(keyword(fr)[c("$P1N", "$P3N")]), chnls)
  
  #update fs
  expect_error(fs[[1]] <- fr, "don't match")
  
  #update colnames for flowSet
  colnames(fs)[c(1,3)] <- chnls
  expect_equal(colnames(fs)[c(1,3)], chnls)
  expect_equivalent(unlist(keyword(fs[[1]])[c("$P1N", "$P3N")]), chnls)
  
  #now [[<- succeeds
  newmarker <- "fsc-h"
  markernames(fr) <- c(A = newmarker)
  fs[[1]] <- fr
  expect_equivalent(markernames(fs[[1]])[1], newmarker)
  
  #swap cols
  colnames(fr)[c(1,3)] <- rev(chnls)
  expect_error(fs[[1]] <- fr, "don't match")
  
})
