library(flowCore)
library(flowUtils)
library(RUnit)
library(XML)

fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InvalidGateTests/FCS_Files","testfile_17parameters.fcs",package="flowCore")

gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InvalidGateTests/Gating-ML_Files","NoDimensionsRectangles.xml",package="flowCore")

source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InvalidGateTests/RUnitScript_Files","performSimpleGateTest.R",package="flowCore"))

test.1 <- function() {
  performSimpleGateTest("CD3",fcsFile,gateFile)
}

test.2 <- function() {
  performSimpleGateTest("RectGate1",fcsFile,gateFile)
}

