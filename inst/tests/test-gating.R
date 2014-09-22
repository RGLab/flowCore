b08 = read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))

#test the construction of rectangle gates using 'rectGate'
test_that("rectGateCreation", {
  expect_is(rectangleGate("FSC-H"=c(125,500)),"rectangleGate")
  expect_is(rectangleGate("FSC-H"=c(125,500),"SSC-H"=c(75,800)),"rectangleGate")
  expect_is(rectangleGate(list("FSC-H"=c(125,500),"SSC-H"=c(75,800))),"rectangleGate")
})

test_that("rectGateUsage", {
	#Test out actually using a gate on some data
    expect_true(isGeneric("%in%"))
	expect_is(b08, "flowFrame")
    
	r = rectangleGate("FSC-H"=c(125,500))
	expect_is(r,"rectangleGate")
	
	expect_true(hasMethod("%in%",c(class(b08),class(r))))
		

	x = getMethod("%in%",c(class(b08),class(r)))(b08,r)
	expect_equal(sum(x), 6906)
})


## test the construction of polygon gates using 'polygonGate'
test_that("polygonGateCreation", {
        vertices <- matrix(c(300,400,600,400,300, 300, 50, 70, 200, 180,150,50), ncol=2)
        colnames(vertices) <- c("FSC-H", "SSC-H")
	    expect_is(polygonGate(.gate=vertices),"polygonGate")
})


## Test gate on some data
test_that("polygonGateUsage", {
	
    vertices <- matrix(c(300,400,600,400,300, 300, 50, 70, 200, 180,150,50), ncol=2)
    colnames(vertices) <- c("FSC-H", "SSC-H")
	p   = polygonGate(vertices)
	expect_is(p,"polygonGate")
	
    expect_true(hasMethod("%in%",c(class(b08),class(p))))
	x   = b08 %in% p
	m = getMethod("%in%",c(class(b08),class(p)))(b08,p)
	expect_equal(sum(x),3785)
    p2   = polygonGate(vertices[-nrow(vertices),])
   	expect_equal(sum(x),sum(b08 %in% p2))
})
