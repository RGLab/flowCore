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

# Load the HowTO.R script in InOutEventsTests folder
source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InOutEventsTests","HowTO.R",package="flowCore"))

# Load the HowTO.R script in InvalidGateTests folder
source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/InvalidGateTests","HowTO.R",package="flowCore"))

# Load the HowTO.R script in ValidGateTests folder
source(system.file("unitTests/standardsCompliance/Gating-ML_ComplianceTests/ValidGateTests","HowTO.R",package="flowCore"))

