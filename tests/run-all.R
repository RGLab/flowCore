library(testthat)
library(flowCore)
#library(methods)
#library(data.table)
#library(utils)

#dataDir <- system.file("extdata",package="flowWorkspaceData")
#resultDir <- system.file("tests/expect_result",package="flowWorkspace")
test_package("flowCore")

expectRes <- readRDS(system.file("tests/expectResults.rds", package = "flowCore"))
#test_file("/home/wjiang2/rglab/workspace/flowWorkspace/examples/InternalTestSuite.R")