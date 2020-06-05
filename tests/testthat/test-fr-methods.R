context("flowFrame methods...")
data("GvHD")
fr <- GvHD[[1]]

test_that("flowFrame constructor", {
  mat <- matrix(1:30,ncol = 3, dimnames = list(NULL, letters[1:3]))
  fr <- flowFrame(mat)
  markernames(fr)[1] <- "A"
  expect_equal(markernames(fr)[1] , c(a = "A"))
})
test_that("range", {
  rng1 <- data.frame("FSC-H" = c(0,1023)
                     ,"SSC-H" = c(0,1023)
                     ,"FL1-H" = c(1,10000)
                     ,"FL2-H" = c(1,10000)
                     ,"FL3-H" = c(1,10000)
                     ,"FL2-A" = c(0,1023)
                     ,"FL4-H" = c(1,10000)
                     ,"Time" = c(0,1023)
                     , row.names = c("min", "max")
                     , check.names = FALSE
                    )
  expect_equal(range(fr), rng1)
  
  expect_equal(range(fr, "instrument"), rng1)
  
  expect_equal(range(fr, type = "instrument"), rng1)
  
  expect_error(range(fr, "FSC-H"), "only accept two")
  
  rng2 <- data.frame("FSC-H" = c(59,1023)
                     ,"SSC-H" = c(6,1023)
                     ,"FL1-H" = c(1,10000)
                     ,"FL2-H" = c(1.000,9221.666)
                     ,"FL3-H" = c(1.000,1131.784)
                     ,"FL2-A" = c(0,1023)
                     ,"FL4-H" = c(1,1162.77)
                     ,"Time" = c(1, 755)
                     , row.names = c("min", "max")
                     , check.names = FALSE
                    )
  expect_equal(range(fr, type = "data")  ,rng2, tolerance = 4e-7)
  expect_equal(range(fr, "data")  ,rng2, tolerance = 4e-7)
  expect_error(range(fr, "FSC-H", type = "data"), "only accept two")
  
  #swap cols of fr
  origcol <- colnames(fr)
  colnames(fr@exprs)[7:8] <- origcol[8:7]
  rg <- as.matrix(range(fr, "data")[7:8])
  rownames(rg) <- NULL
  expect_equal(apply(exprs(fr)[,7:8], 2, range), rg)
})

test_that("fr_remove_redundant_pnx_kw", {
  
expect_equal(length(keyword(fr[,c(1:6,8)])), length(keyword(fr[,-7])))
  
})


test_that("fr_append_cols", {
  
  n <- matrix(1:(nrow(fr)), ncol = 1)
  colnames(n) <- "A"
  m <- matrix(1:(2*nrow(fr)), ncol = 2)
  colnames(m) <- c("B", "C")
  
  # Add single column and make sure min/max keywords set appropriately
  fr_plus <- fr_append_cols(fr, n)
  key_range <- keyword(fr_plus)[c("flowCore_$P9Rmin", "flowCore_$P9Rmax")]
  expect_equal(unname(unlist(key_range)), range(n[,"A"]))
  
  # Add multiple columns
  fr_plus <- fr_append_cols(fr, m)
  key_range <- keyword(fr_plus)[c("flowCore_$P9Rmin", "flowCore_$P9Rmax")]
  expect_equal(unname(unlist(key_range)), range(m[,"B"]))
  key_range <- keyword(fr_plus)[c("flowCore_$P10Rmin", "flowCore_$P10Rmax")]
  expect_equal(unname(unlist(key_range)), range(m[,"C"]))
  
})
