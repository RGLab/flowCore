####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "SimpleBooleanGates.xml" Gating-ML file                                       #
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
# Table of expected results:
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
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","SimpleBooleanGates.xml",package="flowCore")

# The expected results read in from the CSV file
expectedResults <- read.table(csvFile,sep=",",header=TRUE, check.names=FALSE) 

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

test.AndTwoArg0In <- function() {
  performGateTest("AndTwoArg0In", testGates, fcs, expectedResults, numberOfGates)
}

test.AndTwoArg1In <- function() {
  performGateTest("AndTwoArg1In", testGates, fcs, expectedResults, numberOfGates)
}

test.AndTwoArgAllIn <- function() {
  performGateTest("AndTwoArgAllIn", testGates, fcs, expectedResults, numberOfGates)
}

test.AndManyArgs0In <- function() {
  performGateTest("AndManyArgs0In", testGates, fcs, expectedResults, numberOfGates)
}

test.AndManyArgs1In <- function() {
  performGateTest("AndManyArgs1In", testGates, fcs, expectedResults, numberOfGates)
}

test.AndManyArgsSomeIn <- function() {
  performGateTest("AndManyArgsSomeIn", testGates, fcs, expectedResults, numberOfGates)
}

test.AndManyArgsAllIn <- function() {
  performGateTest("AndManyArgsAllIn", testGates, fcs, expectedResults, numberOfGates)
}

test.OrTwoArg0In <- function() {
  performGateTest("OrTwoArg0In", testGates, fcs, expectedResults, numberOfGates)
}

test.OrTwoArg1In <- function() {
  performGateTest("OrTwoArg1In", testGates, fcs, expectedResults, numberOfGates)
}

test.OrTwoArgAllIn <- function() {
  performGateTest("OrTwoArgAllIn", testGates, fcs, expectedResults, numberOfGates)
}

test.OrManyArgs0In <- function() {
  performGateTest("OrManyArgs0In", testGates, fcs, expectedResults, numberOfGates)
}

test.OrManyArgs1In <- function() {
  performGateTest("OrManyArgs1In", testGates, fcs, expectedResults, numberOfGates)
}

test.OrManyArgsSomeIn <- function() {
  performGateTest("OrManyArgsSomeIn", testGates, fcs, expectedResults, numberOfGates)
}

test.OrManyArgsAllIn <- function() {
  performGateTest("OrManyArgsAllIn", testGates, fcs, expectedResults, numberOfGates)
}

test.NotEventInReferencedGate <- function() {
  performGateTest("NotEventInReferencedGate", testGates, fcs, expectedResults, numberOfGates)
}

test.NotEventNotInReferencedGate <- function() {
  performGateTest("NotEventNotInReferencedGate", testGates, fcs, expectedResults, numberOfGates)
}

test.AndMixedTwoArg1In <- function() {
  performGateTest("AndMixedTwoArg1In", testGates, fcs, expectedResults, numberOfGates)
}

test.AndAllDefinedTwoArgAllIn <- function() {
  performGateTest("AndAllDefinedTwoArgAllIn", testGates, fcs, expectedResults, numberOfGates)
}

test.OrMixedTwoArg1In <- function() {
  performGateTest("OrMixedTwoArg1In", testGates, fcs, expectedResults, numberOfGates)
}

test.OrAllDefinedTwoArg0In <- function() {
  performGateTest("OrAllDefinedTwoArg0In", testGates, fcs, expectedResults, numberOfGates)
}

test.NotEventNotInDefinedGate <- function() {
  performGateTest("NotEventNotInDefinedGate", testGates, fcs, expectedResults, numberOfGates)
}

test.XorBothIn <- function() {
  performGateTest("XorBothIn", testGates, fcs, expectedResults, numberOfGates)
}

test.XorBothOut <- function() {
  performGateTest("XorBothOut", testGates, fcs, expectedResults, numberOfGates)
}

test.XorOneIn <- function() {
  performGateTest("XorOneIn", testGates, fcs, expectedResults, numberOfGates)
}

test.XorOtherIn <- function() {
  performGateTest("XorOtherIn", testGates, fcs, expectedResults, numberOfGates)
}

