####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "ParentIdTest.xml" Gating-ML file                                             #
# applied to                                                                       #
#  - "fcs2_int16_50000ev_8par_random.fcs" FCS data file                            #
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
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","ParentIdTest.xml",package="flowCore")

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

test.prt_p01 <- function() {
  performGateTest(gateID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p02 <- function() {
  performGateTest(gateID="prt_p02", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p03 <- function() {
  performGateTest(gateID="prt_p03", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p04_prt_p01 <- function() {
  performGateTest(gateID="prt_p04", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p05_prt_p04 <- function() {
  performGateTest(gateID="prt_p05", parentID="prt_p04", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p06_prt_p01 <- function() {
  performGateTest(gateID="prt_p06", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p07_prt_p06 <- function() {
  performGateTest(gateID="prt_p07", parentID="prt_p06", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p08_prt_p01 <- function() {
  performGateTest(gateID="prt_p08", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_p09_prt_p08 <- function() {
  performGateTest(gateID="prt_p09", parentID="prt_p08", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_01_prt_p01 <- function() {
  performGateTest(gateID="prt_01", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_02_prt_p01 <- function() {
  performGateTest(gateID="prt_02", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_03_prt_p02 <- function() {
  performGateTest(gateID="prt_03", parentID="prt_p02", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_04_prt_p01 <- function() {
  performGateTest(gateID="prt_04", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_05_prt_p01 <- function() {
  performGateTest(gateID="prt_05", parentID="prt_p01", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_06_prt_p03 <- function() {
  performGateTest(gateID="prt_06", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_07_prt_p03 <- function() {
  performGateTest(gateID="prt_07", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_08_prt_p03 <- function() {
  performGateTest(gateID="prt_08", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_09_prt_p03 <- function() {
  performGateTest(gateID="prt_09", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_10_prt_p03 <- function() {
  performGateTest(gateID="prt_10", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_11_prt_p03 <- function() {
  performGateTest(gateID="prt_11", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_12_prt_p03 <- function() {
  performGateTest(gateID="prt_12", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_13_prt_p03 <- function() {
  performGateTest(gateID="prt_13", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_14_prt_p03 <- function() {
  performGateTest(gateID="prt_14", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_15_prt_p03 <- function() {
  performGateTest(gateID="prt_15", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_16_prt_p03 <- function() {
  performGateTest(gateID="prt_16", parentID="prt_p03", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_17_prt_p08 <- function() {
  performGateTest(gateID="prt_17", parentID="prt_p05", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_18_prt_p07 <- function() {
  performGateTest(gateID="prt_18", parentID="prt_p07", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_19_prt_p09 <- function() {
  performGateTest(gateID="prt_19", parentID="prt_p09", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_20_prt_p07 <- function() {
  performGateTest(gateID="prt_20", parentID="prt_p07", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_21_prt_p05 <- function() {
  performGateTest(gateID="prt_21", parentID="prt_p05", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_22_prt_p07 <- function() {
  performGateTest(gateID="prt_22", parentID="prt_p07", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_23_prt_p09 <- function() {
  performGateTest(gateID="prt_23", parentID="prt_p09", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_24_prt_p05 <- function() {
  performGateTest(gateID="prt_24", parentID="prt_p05", testGates, fcs, expectedResults, numberOfGates)
}

test.prt_25_prt_p07 <- function() {
  performGateTest(gateID="prt_25", parentID="prt_p07", testGates, fcs, expectedResults, numberOfGates)
}
