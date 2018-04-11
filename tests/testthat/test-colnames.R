data(GvHD)
fs <- GvHD[1:2]
fr <- GvHD[[1]]

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
  
  #update colnames for flowSet
  colnames(fs)[c(1,3)] <- chnls
  expect_equal(colnames(fs)[c(1,3)], chnls)
  expect_equivalent(unlist(keyword(fs[[1]])[c("$P1N", "$P3N")]), chnls)
})
