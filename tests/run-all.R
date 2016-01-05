library(testthat)
library(flowCore)

expectRes <- readRDS(system.file("tests/expectResults.rds", package = "flowCore"))
test_package("flowCore")


#the FCS files are not included in package due to the size , thus only for internal testing
  #test_file("/home/wjiang2/rglab/workspace/flowCore/inst/tests/IO-testSuite.R")
