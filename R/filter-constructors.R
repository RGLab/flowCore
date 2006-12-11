## last changed: pdh 12-06-06
##    I commended out the last line, which is the uncompleted beginning of a new function

## Happiness is a warm constructor function
## Constructs rectangleGates from the following forms:
##   rectGate("FSC-H"=c(1,100),"SSC-H"=c(100,200))
##   rectGate(list("FSC-H"=c(1,100),"SSC-H"=c(100,200)))
## or 
##   rectGate(m)
## where 'm' is a matrix the form
##     FSC-H SSC-H
## min 1     100
## max 100   200


## ==========================================================================
## rectangleGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rectangleGate <- function(.gate,...,filterId="rectangleGate") {
    if(missing(.gate) || !is.matrix(.gate))
      	.gate <- sapply(if(missing(.gate)) list(...) else .gate,function(x) c("min"=x[1],"max"=x[2]))
	new("rectangleGate",filterId=filterId,parameters=colnames(.gate),min=.gate[1,],max=.gate[2,])
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## ==========================================================================
## polygonGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
polygonGate <- function(filterId="polygonGate", boundaries,...) {
    if(missing(boundaries) || !is.matrix(boundaries))
      boundaries <- sapply(if(missing(boundaries)) list(...) else boundaries, function(x) x)
    new("polygonGate",filterId=filterId, parameters=colnames(boundaries),boundaries=boundaries)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## ==========================================================================
## polytopeGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
polytopeGate <- function(.gate,...,filterId="polytopeGate") {
    if(missing(.gate) || !is.matrix(.gate))
      ##nrowGate <- max(unlist(lapply(list(...),length)))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
         
    new("polytopeGate",filterId=filterId,parameters=colnames(.gate), boundaries=.gate)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## ==========================================================================
## EllipsoideGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ellipsoideGate <- function(.gate, distance,...,filterId="ellipsoidGate") {
    if(missing(.gate) || !is.matrix(.gate))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
      
    new("ellipsoidGate",filterId=filterId,parameters=colnames(.gate), focus=.gate,distance=distance)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## ==========================================================================
## Norm2Filter contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
norm2Filter <- function(x,y,method="covMcd",scale.factor=1,filterId="norm2Gate",...) {
	if(missing(y)) {
		if(length(x)==1)
			stop("You must specify two parameters for a norm2 gate.")
		if(length(x)>2)
			warning("Only the first two parameters will be used.")
		y=x[2]
		x=x[1]
	} else {
		if(length(x)>1 || length(y)>1)
			warning("Only the first two parameters from 'x' and 'y' will be used.")
			x = x[1]
			y = y[1]
	}
	new("norm2Filter",parameters=c(x,y),method=method,scale.factor=scale.factor,filterId=filterId,...)
}
