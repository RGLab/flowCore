context("converting raw vector to integer or numeric that has different bitwidths")
test_that("raw to int", {
  # Integer
  origVec <- c(1:10)
  nPar <- 5
  nByteSize <- 4
  size <- rep(nByteSize, nPar)
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "little")
  # wrong size vec
  expect_error(convertRawBytes(rawVec, isInt = T, colSize = 4, ncol = nPar, isBigEndian = F), "length of 'colSize'")
  
  expect_identical(origVec, convertRawBytes(rawVec, isInt = T, colSize = size, ncol = nPar, isBigEndian = F))
  
  #byte size other than 4
  nByteSize <- 2
  size <- rep(nByteSize, nPar)
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "little")
  expect_identical(origVec, convertRawBytes(rawVec, isInt = T, colSize = size, ncol = nPar, isBigEndian = F))
  
  nByteSize <- 1
  size <- rep(nByteSize, nPar)
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "little")
  expect_identical(origVec, convertRawBytes(rawVec, isInt = T, colSize = size, ncol = nPar, isBigEndian = F))
  
  #big endian
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "big")
  expect_identical(origVec, convertRawBytes(rawVec, isInt = T, colSize = size, ncol = nPar, isBigEndian = T))
  
  #mixed sizes
  size <- c(2,2,2,4,4)
  rawVec1 <- writeBin(origVec[c(1:3,6:8)], raw(), size = 2, endian = "little")
  rawVec2 <- writeBin(origVec[c(4:5,9:10)], raw(), size = 4, endian = "little")
  rawVec <- c(rawVec1[1:6], rawVec2[1:8], rawVec1[7:12], rawVec2[9:16])
  expect_identical(origVec, convertRawBytes(rawVec, isInt = T, colSize = size, ncol = nPar, isBigEndian = F))
  
  
})

test_that("raw to numeric", {
  # double
  origVec <- as.numeric(c(1:10))
  nPar <- 5
  nByteSize <- 8L
  size <- rep(nByteSize, nPar)
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "little")
  expect_identical(origVec,convertRawBytes(rawVec, isInt = F, colSize = size, ncol = nPar, isBigEndian = F))
  
  #big endian
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "big")
  expect_identical(origVec,convertRawBytes(rawVec, isInt = F, colSize = size, ncol = nPar, isBigEndian = T))
  
  nByteSize <- 4L
  size <- rep(nByteSize, nPar)
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "little")
  expect_identical(origVec,convertRawBytes(rawVec, isInt = F, colSize = size, ncol = nPar, isBigEndian = F))
  
  
  
})
