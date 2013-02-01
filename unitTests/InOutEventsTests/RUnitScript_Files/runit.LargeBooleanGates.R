####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "LargeBooleanGates.xml" Gating-ML file                                        #
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
csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-10000_events_random.csv",package="flowCore")

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-10000_events_random.fcs",package="flowCore")

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","LargeBooleanGates.xml",package="flowCore")

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

test.Rect1 <- function() {
  performGateTest("Rect1", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect2 <- function() {
  performGateTest("Rect2", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect3 <- function() {
  performGateTest("Rect3", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect4 <- function() {
  performGateTest("Rect4", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect5 <- function() {
  performGateTest("Rect5", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect6 <- function() {
  performGateTest("Rect6", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect7 <- function() {
  performGateTest("Rect7", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect8 <- function() {
  performGateTest("Rect8", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect9 <- function() {
  performGateTest("Rect9", testGates, fcs, expectedResults, numberOfGates)
}

test.Rect10 <- function() {
  performGateTest("Rect10", testGates, fcs, expectedResults, numberOfGates)
}

test.Bool1 <- function() {
  performGateTest("Bool1", testGates, fcs, expectedResults, numberOfGates)
}
