####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "ComplicatedPolygons.xml" Gating-ML file                                      #
# applied to                                                                       #
#  - "int-15_scatter_events.fcs" FCS data file                                     #
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
csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-15_scatter_events.csv",package="flowCore") 

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-15_scatter_events.fcs",package="flowCore")

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","ComplicatedPolygons.xml",package="flowCore")
 
# The expected results read in from the CSV file
expectedResults <- read.table(csvFile,sep=",",header=TRUE) 

# The FCS data (an instance of the flowFrame class) read in from the FCS data file 
fcs <- read.FCS(fcsFile)

# A list of gates read in from the Gating-ML file
testGates <- read.gatingML(gateFile)

# Change the type of the "fcs@exprs" slot (the actual parameter values) to numeric
fcs@exprs[,1:ncol(fcs@exprs)] <- as.numeric(fcs@exprs)
/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files
# The count of gates within the testgates list
numberOfGates <- xmlSize(testGates)

# Load the performGateTest function
source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/RUnitScript_Files","performGateTest.R",package="flowCore"))

test.G1i1 <- function() {
  performGateTest("G1i1", testGates, fcs, expectedResults, numberOfGates)
}

test.G1i2 <- function() {
  performGateTest("G1i2", testGates, fcs, expectedResults, numberOfGates)
}

test.G1i3 <- function() {
  performGateTest("G1i3", testGates, fcs, expectedResults, numberOfGates)
}

test.G1i4 <- function() {
  performGateTest("G1i4", testGates, fcs, expectedResults, numberOfGates)
}

test.G2i1 <- function() {
  performGateTest("G2i1", testGates, fcs, expectedResults, numberOfGates)
}

test.G2i2 <- function() {
  performGateTest("G2i2", testGates, fcs, expectedResults, numberOfGates)
}

test.G3i1 <- function() {
  performGateTest("G3i1", testGates, fcs, expectedResults, numberOfGates)
}

test.G3i2 <- function() {
  performGateTest("G3i2", testGates, fcs, expectedResults, numberOfGates)
}
test.G3i3 <- function() {
  performGateTest("G3i3", testGates, fcs, expectedResults, numberOfGates)
}

test.G3i4 <- function() {
  performGateTest("G3i4", testGates, fcs, expectedResults, numberOfGates)
}

test.G3i5 <- function() {
  performGateTest("G3i5", testGates, fcs, expectedResults, numberOfGates)
}

test.G3i6 <- function() {
  performGateTest("G3i6", testGates, fcs, expectedResults, numberOfGates)
}

test.G3i7 <- function() {
  performGateTest("G3i7", testGates, fcs, expectedResults, numberOfGates)
}

test.G4i1 <- function() {
  performGateTest("G4i1", testGates, fcs, expectedResults, numberOfGates)
}

test.G4i2 <- function() {
  performGateTest("G4i2", testGates, fcs, expectedResults, numberOfGates)
}

test.G4i3 <- function() {
  performGateTest("G4i3", testGates, fcs, expectedResults, numberOfGates)
}

test.G4i4 <- function() {
  performGateTest("G4i4", testGates, fcs, expectedResults, numberOfGates)
}

test.G5i1 <- function() {
  performGateTest("G5i1", testGates, fcs, expectedResults, numberOfGates)
}

test.G5i2 <- function() {
  performGateTest("G5i2", testGates, fcs, expectedResults, numberOfGates)
}

test.G5i3 <- function() {
  performGateTest("G5i3", testGates, fcs, expectedResults, numberOfGates)
}

test.G5i4 <- function() {
  performGateTest("G5i4", testGates, fcs, expectedResults, numberOfGates)
}

test.G5i5 <- function() {
  performGateTest("G5i5", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i4 <- function() {
  performGateTest("G6i4", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i8 <- function() {
  performGateTest("G6i8", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i16 <- function() {
  performGateTest("G6i16", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i100 <- function() {
  performGateTest("G6i100", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i101 <- function() {
  performGateTest("G6i101", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i1000 <- function() {
  performGateTest("G6i1000", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i1024 <- function() {
  performGateTest("G6i1024", testGates, fcs, expectedResults, numberOfGates)
}

test.G6i2347 <- function() {
  performGateTest("G6i2347", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate6i10000 <- function() {
  performGateTest("G6i10000", testGates, fcs, expectedResults, numberOfGates)
}

test.G7i1 <- function() {
  performGateTest("G7i1", testGates, fcs, expectedResults, numberOfGates)
}
