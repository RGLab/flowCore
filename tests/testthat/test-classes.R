#Make sure our special generics are actually defined as generic.
test_that("generics", {
	expect_true(isGeneric("transform"))
    expect_true(isGeneric("%in%"))
	#Actually resolve a method
    expect_true(hasMethod("%in%",c("flowFrame","rectangleGate")))
    expect_true(hasMethod("%in%",c("flowFrame","rectangleGate")))
    
})