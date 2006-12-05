library(methods)
library(RUnit)
library(flowCore)


#test the construction of rectangle gates using 'rectGate'
test.rectGateCreation = function() {
	checkTrue(is(rectGate("FSC-H"=c(125,500)),"rectangleGate"))
	checkTrue(is(rectGate("FSC-H"=c(125,500),"SSC-H"=c(75,800)),"rectangleGate"))
	checkTrue(is(rectGate(list("FSC-H"=c(125,500),"SSC-H"=c(75,800))),"rectangleGate"))
}

test.rectGateUsage = function() {
	#Test out actually using a gate on some data
	checkTrue(isGeneric("%in%"))
	b08 = read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
	print(class(b08))
	checkTrue(is(b08,"flowFrame"))
	r   = rectGate("FSC-H"=c(125,500))
	checkTrue(is(r,"rectangleGate"))
	print(class(r))
	print(hasMethod("%in%",c(class(b08),class(r))))
	print(getMethod("%in%",c(class(b08),class(r))))	
#	x   = b08 %in% r
	x = getMethod("%in%",c(class(b08),class(r)))(b08,r)
	checkEquals(sum(x),6918)
}