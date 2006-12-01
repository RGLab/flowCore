#test the construction of rectangle gates using 'rectGate'
test.rectGateCreation = function() {
	checkTrue(is(rectGate("FSC-H"=c(125,500)),"rectangleGate"))
	checkTrue(is(rectGate("FSC-H"=c(125,500),"SSC-H"=c(75,800)),"rectangleGate"))
	checkTrue(is(rectGate(list("FSC-H"=c(125,500),"SSC-H"=c(75,800))),"rectangleGate"))
}

test.rectGateUsage = function() {
	require(flowCore)
	#Test out actually using a gate on some data
	b08 = read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
	checkTrue(is(b08,"flowFrame"))
	r   = rectGate("FSC-H"=c(125,500))
	checkTrue(is(r,"rectangleGate"))
	x   = b08 %in% r
	checkEquals(sum(x),6918)
}