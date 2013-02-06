####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#                                                                                  #
# File includes the performGateTest function                                       #
# Parametes:                                                                       #
# - gateID - the ID of the gate; corresponds to filterId slot the particular gate  #
#            object and to the column name of the table with expected results      #
# - testGates - list of gates                                                      #
# - fcs - the FCS data as an instance of the flowFrame                             #
# - expectedResults - table of expected results                                    #
#   - columns correspond to events in FCS file                                     #
#   - rows correspond to events in FCS file                                        #
#   - individual cells: 1 = the event is expected to be found in the gate          #
#                       0 = the event is NOT expected to be found in the gate      #
# - numberOfGates - the number of gates in testgates object, default is 0, which   #
#                   means that the functions computes the real numberOfGates.      #
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
# (list of unknown gates) and is applied on the FCS data. The results are compared against
# the appropriate column of the expectedResults table. Any discrepancy between the computed
# and expected results (gate membership of each event) is reported within the RUnit framework.
#
# Parametes:
# - gateID - the ID of the gate; corresponds to filterId slot the particular gate object
#            and to the column name of the table with expected results
# - testGates - list of gates
# - fcs - the FCS data as an instance of the flowFrame
# - expectedResults - table of expected results
#   - columns correspond to events in FCS file
#   - rows correspond to events in FCS file
#   - individual cells: 1 = the event is expected to be found in the gate
#                       0 = the event is NOT expected to be found in the gate
# - numberOfGates - the number of gates in testgates object, default is 0, which means
#                   that the functions computes the real numberOfGates.
performGateTest <- function(gateID, testGates, fcs, expectedResults, numberOfGates=0) {
   # Indicated if we have found the gateID in the list of all gates,
   # start as FALSE and switch to TRUE once the gate has been found
   gateFound <- FALSE

   # numberOfGates == 0 means that it was not provided and we compute the number here
   if(numberOfGates == 0) numberOfGates <- xmlSize(testGates)

   # Search all gates for gateID and proceed with the body if gate was found
   for(gateIndex in 1:numberOfGates) if(testGates[[gateIndex]]@filterId==gateID)
   {
      # Calculate the actual 0/1 vector indicating for each event whether it is in the gate
      detected <- as.numeric(getMethod("%in%",c(class(fcs),class(testGates[[gateIndex]])))(fcs,testGates[[gateIndex]]))
      # Also remember that the gate has been found
      gateFound <- TRUE
   }
   
   # Check if we have found the gate and fail this test if not; also, report the appropriate error message
   message <- paste("Gate", gateID, "not found.")
   checkTrue(gateFound, message)

   # Get the column of expected result correspond to gateID
   expected <- expectedResults[,gateID]
   length1 <- length(expected)
   length2 <- length(detected)
   
   # Verify that the number of events is the same in both FCS and expected result file;
   # Also, if not then fail this test with the appropriate error message.
   checkTrue((length1 == length2), paste("Number of events in FCS does not correspond to number of events in expected results file, for",gateID,sep=""))

   # Prepare a message just in case we detect a discrepancy for a particular event
   message <- "Wrong event numbers:"
   
   # No discrepancy has been detected so far, we switch this once we detect one
   eventsCorrect <- TRUE
   
   # Go through all events and compare the detected vs. expected value for each event
   for(eventIndex in 1:length1)
   {
      if(expected[eventIndex] != detected[eventIndex]) {
        # In case there is a discrepancy, add the event number to the message
        # and remember that there is at least one error
        message <- paste(message, eventIndex)
        eventsCorrect <- FALSE
      }
   }

   # Fail the test if there was at least one discrepancy between the
   # expected and detected results; report the "wrong" event numbers in this case.
   checkTrue(eventsCorrect, message)
}
