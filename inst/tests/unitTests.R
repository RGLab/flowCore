## A simple RUnit test harness
package = "flowCore"

if(require("RUnit",quietly=TRUE)) {
	path = if(file.exists(file.path("..","unitTests"))) file.path("..","unitTests") else system.file("unitTests",package=package)
	if(require(package,character.only=TRUE,quietly=TRUE)) {
		testSuite = defineTestSuite(name=paste(package,"unit testting"),dirs=path)
		tests = runTestSuite(testSuite)
		printTextProtocol(tests)
#		if(tests$nFail > 0)
#			stop("Some tests failed to complete.")
	}
}