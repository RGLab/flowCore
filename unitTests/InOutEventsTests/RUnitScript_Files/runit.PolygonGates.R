####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "PolygonGates.xml" Gating-ML file   runTestFile(PolygonGates)                 # 
# applied to                                                                       #
#  - "int-gating_test_file.fcs" FCS data file                                      #
#                                                                                  #
# @author Stanley Wong, Josef Spidlen, Ryan Brinkman                               #
#                                                                                  #
# File distributed under terms of the GNU Lesser General Public License (LGPL)     #
####################################################################################

# Ensure loading of the required libraries
library(flowCore)
library(flowUtils)
library(RUnit)
library(XML)

# Variables shared for all the functions within this script file
#
# Table of expected results:runTestFile(PolygonGates)
# - columns correspond to gate IDs
# - rows correspond to events in FCS file
# - individual cells: 1 = the event is expected to be found in the gate
#                     0 = the event is NOT expected to be found in the gate
#####################################################################################################################
# Sinlge Event

# csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-gating_test_file.csv",package="flowCore")

# Full file name of the FCS data file used in tests within this script file
# fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-gating_test_file.fcs",package="flowCore")
#####################################################################################################################  
# Two Events

csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-gating_test_file_2ev.csv",package="flowCore") 

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-gating_test_file_2ev.fcs",package="flowCore")
#####################################################################################################################

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","PolygonGates.xml",package="flowCore")

# The expected results read in from the CSV file
expectedResults <- read.table(csvFile,sep=",",header=TRUE) 

# The FCS data (an instance of the flowFrame class) read in from the FCS data file 
fcs <- read.FCS(fcsFile)

# A list of gates read in from the Gating-ML file
testGates <- read.gatingML(gateFile)

# Change the type of the "fcs@exprs" slot (the actual parameter values) to numeric
fcs@exprs[,1:ncol(fcs@exprs)] <- as.numeric(fcs@exprs)

# The count of gates within the testgates list
numberOfGates <- xmlSize(testGates)

# Load the performGateTest function
source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/RUnitScript_Files","performGateTest.R",package="flowCore"))

test.RectangleInside <- function() {
  performGateTest("RectangleInside", testGates, fcs, expectedResults, numberOfGates)
}

test.RectangleBoundary <- function() {
  performGateTest("RectangleBoundary", testGates, fcs, expectedResults, numberOfGates)
}

test.RectangleOutside <- function() {
  performGateTest("RectangleOutside", testGates, fcs, expectedResults, numberOfGates)
}

test.SimpleConcaveInside <- function() {
  performGateTest("SimpleConcaveInside", testGates, fcs, expectedResults, numberOfGates)
}

test.SimpleConcaveBoundary <- function() {
  performGateTest("SimpleConcaveBoundary", testGates, fcs, expectedResults, numberOfGates)
}

test.SimpleConcaveOutside <- function() {
  performGateTest("SimpleConcaveOutside", testGates, fcs, expectedResults, numberOfGates)
}

test.ConcaveInside <- function() {
  performGateTest("ConcaveInside", testGates, fcs, expectedResults, numberOfGates)
}

test.ConcaveBoundary <- function() {
  performGateTest("ConcaveBoundary", testGates, fcs, expectedResults, numberOfGates)
}

test.ConcaveOutside <- function() {
  performGateTest("ConcaveOutside", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleTopInside <- function() {
  performGateTest("NonSimpleTopInside", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleBottomInside <- function() {
  performGateTest("NonSimpleBottomInside", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleLeftOutside <- function() {
  performGateTest("NonSimpleLeftOutside", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleRightOutside <- function() {
  performGateTest("NonSimpleRightOutside", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleBoundaryCrossingPoint <- function() {
  performGateTest("NonSimpleBoundaryCrossingPoint", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleBoundaryForwardSlant <- function() {
  performGateTest("NonSimpleBoundaryForwardSlant", testGates, fcs, expectedResults, numberOfGates)
}

test.NonSimpleBoundaryBackSlant <- function() {
  performGateTest("NonSimpleBoundaryBackSlant", testGates, fcs, expectedResults, numberOfGates)
}

test.ComplicatedNonSimple1 <- function() {
  performGateTest("ComplicatedNonSimple1", testGates, fcs, expectedResults, numberOfGates)
}

test.ComplicatedNonSimple2 <- function() {
  performGateTest("ComplicatedNonSimple2", testGates, fcs, expectedResults, numberOfGates)
}

test.ComplicatedNonSimple3 <- function() {
  performGateTest("ComplicatedNonSimple3", testGates, fcs, expectedResults, numberOfGates)
}

test.ComplicatedNonSimple4 <- function() {
  performGateTest("ComplicatedNonSimple4", testGates, fcs, expectedResults, numberOfGates)
}








