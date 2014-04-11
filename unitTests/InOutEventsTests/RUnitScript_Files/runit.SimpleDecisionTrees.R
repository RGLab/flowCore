####################################################################################
# This R script file is part RUnit tsets of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# Tests of expected vs. detected events in gates of                                #
#  - "SimpleDecisionTrees.xml" Gating-ML file                                      #
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
gateFile <- system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests/Gating-ML_Files","SimpleDecisionTrees.xml",package="flowCore")

# The expected results read in from the CSV file
expectedResults <- read.table(csvFile,sep=",",header=TRUE,check.names=FALSE) 

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

test.G1D-L <- function() {
  performGateTest("G1D-L")
}

test.G1D-G <- function() {
  performGateTest("G1D-G")
}

test.G1D-E <- function() {
  performGateTest("G1D-E")
}

test.G2D-LL <- function() {
  performGateTest("G2D-LL")
}

test.G2D-LG <- function() {
  performGateTest("G2D-LG")
}

test.G2D-GL <- function() {
  performGateTest("G2D-GL")
}

test.G2D-GG <- function() {
  performGateTest("G2D-GG")
}

test.G4D-LLLL <- function() {
  performGateTest("G4D-LLLL")
}

test.G4D-LLLG <- function() {
  performGateTest("G4D-LLLG")
}

test.G4D-LLGL <- function() {
  performGateTest("G4D-LLGL")
}

test.G4D-LLGG <- function() {
  performGateTest("G4D-LLGG")
}

test.G4D-LGLL <- function() {
  performGateTest("G4D-LGLL")
}

test.G4D-LGLG <- function() {
  performGateTest("G4D-LGLG")
}

test.G4D-LGGL <- function() {
  performGateTest("G4D-LGGL")
}

test.G4D-LGGG <- function() {
  performGateTest("G4D-LGGG")
}

test.G4D-GLLL <- function() {
  performGateTest("G4D-GLLL")
}

test.G4D-GLLG <- function() {
  performGateTest("G4D-GLLG")
}

test.G4D-GLGL <- function() {
  performGateTest("G4D-GLGL")
}

test.G4D-GLGG <- function() {
  performGateTest("G4D-GLGG")
}

test.G4D-GGLL <- function() {
  performGateTest("G4D-GGLL")
}

test.G4D-GGLG <- function() {
  performGateTest("G4D-GGLG")
}

test.G4D-GGGL <- function() {
  performGateTest("G4D-GGGL")
}

test.G4D-GGGG <- function() {
  performGateTest("G4D-GGGG")
}

