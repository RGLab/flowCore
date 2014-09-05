library(testthat)
library(flowCore)
#library(methods)
#library(data.table)
#library(utils)

#dataDir <- system.file("extdata",package="flowWorkspaceData")
#resultDir <- system.file("tests/expect_result",package="flowWorkspace")
test_package("flowCore")

expectRes <- readRDS(system.file("tests/expectResults.rds", package = "flowCore"))
#the FCS files are not included in package due to the size , thus only for internal testing
#test_file(system.file("tests/IO-testSuite.R", package = "flowCore"))
