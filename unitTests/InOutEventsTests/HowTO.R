####################################################################################
# This R script file is part RUnit tests of Bioconductor/flowCore package intended #
# to test Gating-ML compatibility.                                                 #
#										   # 
# Test of all the functions defined by the test suite                              #
#  - the expression for the test files                                             # 
#    - start with runit.                                                           #
#    - end with .r or .R                                                           #
#    - example: runit.xxx.r or runit.xxx.R                                         #
#  - the expression for the test functions                                         #
#    - start with test.                                                            #
#    - example: test.xxx                                                           #
#                                                                                  #
# @author Stanley Wong, Josef Spidlen, Ryan Brinkman                               #
#                                                                                  #
# File distributed under terms of the GNU Lesser General Public License (LGPL)     #
####################################################################################

# Ensures loading of the required libraries
library(flowCore)
library(flowUtils)
library(XML)
library(RUnit)

# Defines a test suite
testsuite <- defineTestSuite(
  "GatingTestSuite",
  dir=system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests","RUnitScript_Files",package="flowCore"),
  testFileRegexp="^runit.+\\.[rR]$", 
  testFuncRegexp="^test.+")

# Runs all the test functions defined by the test suite
testResult <- runTestSuite(testsuite)

# Prints a protocol of a test run
#summary(testResult)

printHTMLProtocol(testResult,file="InOutResult.html")