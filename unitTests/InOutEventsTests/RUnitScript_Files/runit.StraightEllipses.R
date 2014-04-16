####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "LargeEllipsoids.xml" Gating-ML file                                          #
# applied to                                                                       #
#  - "int-10000_events_random.fcs" FCS data file                                   #
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
csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","fcs2_int16_50000ev_8par_random.csv",package="flowCore")

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","fcs2_int16_50000ev_8par_random.fcs",package="flowCore")

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","StraightEllipses.xml",package="flowCore")


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

test.StraightEllipse_1 <- function() {
  performGateTest("StraightEllipse_1", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_2 <- function() {
  performGateTest("StraightEllipse_2", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_3 <- function() {
  performGateTest("StraightEllipse_3", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_4 <- function() {
  performGateTest("StraightEllipse_4", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_5 <- function() {
  performGateTest("StraightEllipse_5", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_6 <- function() {
  performGateTest("StraightEllipse_6", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_7 <- function() {
  performGateTest("StraightEllipse_7", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_8 <- function() {
  performGateTest("StraightEllipse_8", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_9 <- function() {
  performGateTest("StraightEllipse_9", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_10 <- function() {
  performGateTest("StraightEllipse_10", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_11 <- function() {
  performGateTest("StraightEllipse_11", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_12 <- function() {
  performGateTest("StraightEllipse_12", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_13 <- function() {
  performGateTest("StraightEllipse_13", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_14 <- function() {
  performGateTest("StraightEllipse_14", testGates, fcs, expectedResults, numberOfGates)
}

test.StraightEllipse_15 <- function() {
  performGateTest("StraightEllipse_15", testGates, fcs, expectedResults, numberOfGates)
}



