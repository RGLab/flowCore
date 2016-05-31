context("sort bytes that are stored as mixed endian")
origVec <- c(1:10)
origVec_num <- as.numeric(c(1:10))

test_that("raw to int", {
  # Integer
  
  
  nByteSize <- 4
  rawVec <- writeBin(origVec, raw(), size = nByteSize, endian = "little")
  byte_order <- c(3,2,0,1)

  all.equal(sortBytes(rawVec, byte_order), sortBytes1(rawVec, byte_order))
})
