####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "LargeRectangleTest.xml" Gating-ML file                                       #
# applied to                                                                       #
#  - "int-homogenous_matrix.fcs" FCS data file                                     #
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
csvFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/ExpectedResult_Files","int-homogenous_matrix.csv",package="flowCore") 

# Full file name of the FCS data file used in tests within this script file
fcsFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/FCS_Files","int-homogenous_matrix.fcs",package="flowCore")

# Full file name of the Gating-ML file used in tests within this script file
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","LargeRectangleTest.xml",package="flowCore")

# The expected results read in from the CSV file
expectedResults <- read.table(csvFile,sep=",",header=TRUE) 

# The FCS data (an instance of the flowFrame class) read in from the FCS data file 
fcs <- read.FCS(fcsFile)

# A list of gates read in from the Gating-ML file
testGates <- read.gatingML(gateFile)

# Change the type of the "fcs@exprs" slot (the actual parame/Gating-ML_ComplianceTests/InOutEventsTests/RUnitScript_Filester values) to numeric
fcs@exprs[,1:ncol(fcs@exprs)] <- as.numeric(fcs@exprs)

# The count of gates within the testgates list
numberOfGates <- xmlSize(testGates)

# Load the performGateTest function
source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/RUnitScript_Files","performGateTest.R",package="flowCore"))

test.Gate01 <- function() {
  performGateTest("Gate01", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate02 <- function() {
  performGateTest("Gate02", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate03 <- function() {
  performGateTest("Gate03", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate04 <- function() {
  performGateTest("Gate04", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate05 <- function() {
  performGateTest("Gate05", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate06 <- function() {
  performGateTest("Gate06", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate07 <- function() {
  performGateTest("Gate07", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate08 <- function() {
  performGateTest("Gate08", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate11 <- function() {
  performGateTest("Gate11", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate12 <- function() {
  performGateTest("Gate12", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate13 <- function() {
  performGateTest("Gate13", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate14 <- function() {
  performGateTest("Gate14", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate15 <- function() {
  performGateTest("Gate15", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate16 <- function() {
  performGateTest("Gate16", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate17 <- function() {
  performGateTest("Gate17", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate18 <- function() {
  performGateTest("Gate18", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate21 <- function() {
  performGateTest("Gate21", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate22 <- function() {
  performGateTest("Gate22", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate23 <- function() {
  performGateTest("Gate23", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate24 <- function() {
  performGateTest("Gate24", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate25 <- function() {
  performGateTest("Gate25", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate26 <- function() {
  performGateTest("Gate26", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate27 <- function() {
  performGateTest("Gate27", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate28 <- function() {
  performGateTest("Gate28", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate29 <- function() {
  performGateTest("Gate29", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate30 <- function() {
  performGateTest("Gate30", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate31 <- function() {
  performGateTest("Gate31", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate32 <- function() {
  performGateTest("Gate32", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate33 <- function() {
  performGateTest("Gate33", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate34 <- function() {
  performGateTest("Gate34", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate35 <- function() {
  performGateTest("Gate35", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate36 <- function() {
  performGateTest("Gate36", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate37 <- function() {
  performGateTest("Gate37", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate38 <- function() {
  performGateTest("Gate38", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate39 <- function() {
  performGateTest("Gate39", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate40 <- function() {
  performGateTest("Gate40", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate41 <- function() {
  performGateTest("Gate41", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate42 <- function() {
  performGateTest("Gate42", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate43 <- function() {
  performGateTest("Gate43", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate44 <- function() {
  performGateTest("Gate44", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate45 <- function() {
  performGateTest("Gate45", testGates, fcs, expectedResults, numberOfGates)
}

test.Gate46 <- function() {
  performGateTest("Gate46", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep1 <- function() {
  performGateTest("Gatep1", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep2 <- function() {
  performGateTest("Gatep2", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep3 <- function() {
  performGateTest("Gatep3", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep4 <- function() {
  performGateTest("Gatep4", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep5 <- function() {
  performGateTest("Gatep5", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep6 <- function() {
  performGateTest("Gatep6", testGates, fcs, expectedResults, numberOfGates)
}

test.Gatep7 <- function() {
  performGateTest("Gatep7", testGates, fcs, expectedResults, numberOfGates)
}






































