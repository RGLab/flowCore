#Make sure our special generics are actually defined as generic.
test.generics = function() {
	checkTrue(isGeneric("transform"))
	checkTrue(isGeneric("%in%"))
	#Actually resolve a method
	checkTrue(hasMethod("%in%",c("flowFrame","rectangleGate")))
	print(hasMethod("%in%",c("flowFrame","rectangleGate")))
	print(getMethods("%in%"))
	print(get("%in%"))
}