####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# File includes the performSimpleGateTest and testFun functions                    #
# Parametes:                                                                       #
# - gateID - the ID of the gate; corresponds to filterId slot the particular gate  #
#            object                                                                #
# - fcsFile - the FCS data file used in tests                                      #
# - gateFile - the Gating-ML file used in tests                                    #
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

# Function performing the actual test of a single gate.
#
# The gate to be tested is identified by gateID. The gate is taken from the testgates list
# (list of unknown gates) and is applied on the FCS data. Here the checkException function
# is used to test if the Gating-ML function works properly. The result is reported within
# the RUnit framework.
#
# Parametes:
# - gateID - the ID of the gate; corresponds to filterId slot the particular gate object
# - fcsFile - the FCS data file used in tests
# - gateFile - the Gating-ML file used in tests
performSimpleGateTest <- function(gateID, fcsFile, gateFile)
{
   message <- paste("read.gatingML doesn't work!!!")
   checkException(testFun(gateID, fcsFile, gateFile),message,silent=TRUE)
}

testFun <- function(gateID,fcsFile,gateFile) {
   
   # The FCS data (an instance of the flowFrame class) read in from the FCS data file 
   fcs <- read.FCS(fcsFile)
   
   # Check if fcs is actual a FCS file
   message <- paste("Not a FCS file!!!")
   FCSFile <- is(fcs,"flowFrame")
   checkTrue(FCSFile,message)
   
   # Change the type of the "fcs@exprs" slot (the actual parameter values) to numeric
   fcs@exprs[,1:ncol(fcs@exprs)] <- as.numeric(fcs@exprs)
   
   # A list of gates read in from the Gating-ML file
   testGates <- read.gatingML(gateFile)

   # The count of gates within the testgates list
   numberOfGates <- xmlSize(testGates)
   
   # Indicated if we have found the gateID in the list of all gates,
   # start as FALSE and switch to TRUE once the gate has been found
   gateFound <- FALSE
   
   # Search all gates for gateID and proceed with the body if gate was found
   for(gateIndex in 1:numberOfGates) if(testGates[[gateIndex]]@filterId==gateID)
   {
      # Calculate the actual 0/1 vector indicating for each event whether it is in the gate
      detected <- as.numeric(getMethod("%in%",c(class(fcs),class(testGates[[gateIndex]])))(fcs,testGates[[gateIndex]]))
      # Check if getMethod works properly
      message <- paste("getMethod doesn't work!!!")
      checkTrue(is(detected,"numeric"))
      # Also remember that the gate has been found
      gateFound <- TRUE
   }

   # Check if we have found the gate and fail this test if not; also, report the appropriate error message
   message <- paste("Gate", gateID, "not found.")
   checkTrue(gateFound, message)
   
   # Get the total number of events
   length = length(detected)
   message <- "Events inside the gate:"
 
   for(eventIndex in 1:length) if(detected[eventIndex]==1)
   {
      message <- paste(message, eventIndex)
   }
   cat(message,"\n")
}



