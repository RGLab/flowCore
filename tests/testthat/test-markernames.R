data(GvHD)
fs <- GvHD[1:2]
fr <- GvHD[[1]]

test_that("markernames", {
  markers <- c('CD15 FITC','CD45 PE','CD14 PerCP','CD33 APC')
  expect_equal(markernames(fr), markers)
  
  expect_equal(markernames(fs), markers)
  
  markers.new <- c("15", "14")
  names(markers.new) <- c("FL1-H", "FL3-H")
  markernames(fs[[1]]) <- markers.new
  
  expect_warning(res <- markernames(fs), "not consistent")
  expect_equal(unique(fsApply(fs, markernames, simplify = FALSE)), res)
})

test_that("markernames<-", {
  chnls <- c("FL1-H", "FL3-H")
  markers <- c("CD15", "CD14")
  names(markers) <- chnls
  #invalid type
  expect_error(markernames(fr) <- data.frame(name = chnls, desc = markers) , "named character")
  
  #update markers for flowFrame
  markernames(fr) <- markers
  expect_equivalent(markernames(fr)[c(1,3)], markers)
  expect_equivalent(unlist(keyword(fr)[c("$P3S", "$P5S")]), markers)
  
  #update markers for flowSet
  markernames(fs) <-  markers
  expect_equivalent(markernames(fs)[c(1,3)], markers)
  expect_equivalent(unlist(keyword(fr)[c("$P3S", "$P5S")]), markers)
  
  #incorrect channel name
  names(markers)[1] <- "FL1"
  expect_error(markernames(fr) <- markers , "not found")
  
  #NA value in names
  names(markers)[1] <- NA
  expect_error(markernames(fr) <- markers , "not found")
  
  #not named character
  names(markers) <- NULL
  expect_error(markernames(fr) <- markers , "named character")
  
  
})
