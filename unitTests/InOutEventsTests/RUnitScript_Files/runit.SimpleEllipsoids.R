####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "SimpleEllipsoids.xml" Gating-ML file                                         #
# applied to                                                                       #
#  - "int-gating_test_file_4D.fcs" FCS data file                                   #
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

# csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-gating_test_file_4D.csv",package="flowCore") 

# Full file name of the FCS data file used in tests within this script file
# fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-gating_test_file_4D.fcs",package="flowCore")
#####################################################################################################################  
# Two Events

csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-gating_test_file_4D_2ev.csv",package="flowCore") 

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-gating_test_file_4D_2ev.fcs",package="flowCore")
#####################################################################################################################

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","SimpleEllipsoids.xml",package="flowCore")

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

test.WithinDistanceLimit <- function() {
  performGateTest("WithinDistanceLimit", testGates, fcs, expectedResults, numberOfGates)
}

test.WithinDistanceLimitSmall <- function() {
  performGateTest("WithinDistanceLimitSmall", testGates, fcs, expectedResults, numberOfGates)
}

test.EqualDistanceLimit <- function() {
  performGateTest("EqualDistanceLimit", testGates, fcs, expectedResults, numberOfGates)
}

test.OutsideDistanceLimitSmall <- function() {
  performGateTest("OutsideDistanceLimitSmall", testGates, fcs, expectedResults, numberOfGates)
}

test.OutsideDistanceLimit <- function() {
  performGateTest("OutsideDistanceLimit", testGates, fcs, expectedResults, numberOfGates)
}

test.CloserFoci <- function() {
  performGateTest("CloserFoci", testGates, fcs, expectedResults, numberOfGates)
}

test.CloserFociSmall <- function() {
  performGateTest("CloserFociSmall", testGates, fcs, expectedResults, numberOfGates)
}

test.FartherFociSmall <- function() {
  performGateTest("FartherFociSmall", testGates, fcs, expectedResults, numberOfGates)
}

test.FartherFoci <- function() {
  performGateTest("FartherFoci", testGates, fcs, expectedResults, numberOfGates)
}

test.WithinDistanceLimit4D <- function() {
  performGateTest("WithinDistanceLimit4D", testGates, fcs, expectedResults, numberOfGates)
}

test.WithinDistanceLimit4DSmall <- function() {
  performGateTest("WithinDistanceLimit4DSmall", testGates, fcs, expectedResults, numberOfGates)
}

test.EqualDistanceLimit4D <- function() {
  performGateTest("EqualDistanceLimit4D", testGates, fcs, expectedResults, numberOfGates)
}

test.OutsideDistanceLimit4DSmall <- function() {
  performGateTest("OutsideDistanceLimit4DSmall", testGates, fcs, expectedResults, numberOfGates)
}

test.OutsideDistanceLimit4D <- function() {
  performGateTest("OutsideDistanceLimit4D", testGates, fcs, expectedResults, numberOfGates)
}

test.CloserFoci4D <- function() {
  performGateTest("CloserFoci4D", testGates, fcs, expectedResults, numberOfGates)
}

test.CloserFociSmall4D <- function() {
  performGateTest("CloserFociSmall4D", testGates, fcs, expectedResults, numberOfGates)
}

test.FartherFociSmall4D <- function() {
  performGateTest("FartherFociSmall4D", testGates, fcs, expectedResults, numberOfGates)
}

test.FartherFoci4D <- function() {
  performGateTest("FartherFoci4D", testGates, fcs, expectedResults, numberOfGates)
}

