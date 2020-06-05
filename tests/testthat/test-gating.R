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
	expect_equal(sum(x), 6922)
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
   	
   	#create edge cases
   	vertices <- matrix(c(300,300
   	                     ,600,300
   	                     ,800, 100
   	                     , 100, 100
   	                     ), byrow = TRUE,ncol=2)
   	colnames(vertices) <- c("FSC-H", "SSC-H")
   	p   = polygonGate(vertices)

   	#simulate data by adding them to the polygon edges   	
   	fr <- flowFrame(vertices)
   	x   = fr %in% p
   	expect_equal(sum(x),4)
   	
   	#mv one in and one out
   	fr@exprs[1,] <- fr@exprs[1,] + 100
   	fr@exprs[2,] <- fr@exprs[2,] - 100
   	x   = fr %in% p
   	expect_equal(sum(x), 3)
   	
   	
   	#more along the edges
   	fr@exprs[1,2] <- fr@exprs[1,2] - 100
   	fr@exprs[2,2] <- fr@exprs[2,2] + 100
   	x   = fr %in% p
   	expect_equal(sum(x), 4)
   	
   	

   	#add two more
   	pts <- fr@exprs[1,]
   	pts[2] <- pts[2] - 20
   	fr@exprs <- rbind(fr@exprs, pts)
   	pts[2] <- pts[2] - 20
   	fr@exprs <- rbind(fr@exprs, pts)
   	x   = fr %in% p
   	expect_equal(sum(x), 6)
   	
   	#another gate
   	vertices <- matrix(c(300,300
   	                     ,600,300
   	                     ,800, 200
   	                     ,650,100
   	                     ,310, 100
   	                     , 100, 200
   	), byrow = TRUE,ncol=2)
   	colnames(vertices) <- c("FSC-H", "SSC-H")
   	p   = polygonGate(vertices)
   	fr <- flowFrame(vertices)
   	pts <- fr@exprs[1,]
   	pts[2] <- pts[2] - 20
   	fr@exprs <- rbind(fr@exprs, pts)
   	pts <- fr@exprs[2,]
   	pts[2] <- pts[2] - 20
   	fr@exprs <- rbind(fr@exprs, pts)
   	x   = fr %in% p
   	expect_equal(sum(x),8)
   	
   	#add in-cell event that happens to intersect with the right vertex
   	pts <- fr@exprs[2,]
   	pts[2] <- vertices[3,2]
   	fr@exprs <- rbind(fr@exprs, pts)
   	x   = fr %in% p
   	expect_equal(sum(x),9)
   	
   	#add out-cell event that happens to intersect with the right vertex
   	pts[1] <- vertices[6,1] - 50
   	fr@exprs <- rbind(fr@exprs, pts)
   	p@boundaries[6,2] <- p@boundaries[6,2] + 50
   	x   = fr %in% p
   	expect_equal(sum(x),8)
   	
   	# plot(rbind(p@boundaries, fr@exprs), type = "n")
   	# polygon(p@boundaries)
   	# points(fr@exprs[x,], col = "red")
   	# points(fr@exprs[!x,,drop = F], col = "blue")
   	# title(sum(x))
})
#polytopeGate will be deprecated because even gatingML2 no longer supports it
