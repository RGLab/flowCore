####################################################################################
# This R script file is part RUnit tsets of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "PolytopeGates.xml" Gating-ML file                                         #
# applied to                                                                       #
#  - "int-gating_test_file_4D.fcs" FCS data file                                      #
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
csvFile <- system.file("unitTests/gating","int-gating_test_file_4D_2ev.csv",package="flowCore") 

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/gating","int-gating_test_file_4D_2ev.fcs",package="flowCore")

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/gating","PolytopeGates.xml",package="flowCore")

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
source(system.file("unitTests/gating","performGateTest.R",package="flowCore"))

test.LessThanMin <- function() {
  performGateTest("LessThanMin", testGates, fcs, expectedResults, numberOfGates)
}

test.EqualToMin <- function() {
  performGateTest("EqualToMin", testGates, fcs, expectedResults, numberOfGates)
}

test.GreaterThanMin <- function() {
  performGateTest("GreaterThanMin", testGates, fcs, expectedResults, numberOfGates)
}

test.LessThanMax <- function() {
  performGateTest("LessThanMax", testGates, fcs, expectedResults, numberOfGates)
}

test.EqualToMax <- function() {
  performGateTest("EqualToMax", testGates, fcs, expectedResults, numberOfGates)
}

test.GreaterThanMax <- function() {
  performGateTest("GreaterThanMax", testGates, fcs, expectedResults, numberOfGates)
}

test.BetweenMinAndMax <- function() {
  performGateTest("BetweenMinAndMax", testGates, fcs, expectedResults, numberOfGates)
}

test.EqualToMinAndLessThanMax <- function() {
  performGateTest("EqualToMinAndLessThanMax", testGates, fcs, expectedResults, numberOfGates)
}

test.GreaterThanMinAndEqualToMax <- function() {
  performGateTest("GreaterThanMinAndEqualToMax", testGates, fcs, expectedResults, numberOfGates)
}

test.EqualToMinAndMax <- function() {
  performGateTest("EqualToMinAndMax", testGates, fcs, expectedResults, numberOfGates)
}

test.GreaterThanMinAndMax <- function() {
  performGateTest("GreaterThanMinAndMax", testGates, fcs, expectedResults, numberOfGates)
}

test.LessThanMinAndMax <- function() {
  performGateTest("LessThanMinAndMax", testGates, fcs, expectedResults, numberOfGates)
}

test.MinGreaterThanMax <- function() {
  performGateTest("MinGreaterThanMax", testGates, fcs, expectedResults, numberOfGates)
}

test.InNoDimensions <- function() {
  performGateTest("InNoDimensions", testGates, fcs, expectedResults, numberOfGates)
}

test.InOneDimensions <- function() {
  performGateTest("InOneDimensions", testGates, fcs, expectedResults, numberOfGates)
}

test.InAllDimensions <- function() {
  performGateTest("InAllDimensions", testGates, fcs, expectedResults, numberOfGates)
}




