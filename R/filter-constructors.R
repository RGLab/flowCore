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
rectangleGate <- function(filterId="rectangleGate",.gate,...) {
    if(missing(.gate) || !is.matrix(.gate))
	.gate <- sapply(if(missing(.gate)) list(...) else .gate,function(x) c("min"=x[1],"max"=x[2]))
	new("rectangleGate",filterId=filterId,parameters=colnames(.gate),min=.gate[1,],max=.gate[2,])
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## ==========================================================================
## polygonGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
polygonGate <- function(filterId="polygonGate",.gate,...) {
	if(missing(.gate) || !is.matrix(.gate))
	.gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) c("min"=x[1],"max"=x[2]))

        new("polygonGate",filterId=filterId,parameters=colnames(.gate),boundaries=.gate)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## ==========================================================================
## polytopeGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
polytopeGate <- function(filterId="polytopeGate",.gate,...) {
    if(missing(.gate) || !is.matrix(.gate))
      ##nrowGate <- max(unlist(lapply(list(...),length)))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
         
    new("polytopeGate",filterId=filterId,parameters=colnames(.gate), boundaries=.gate)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## ==========================================================================
## EllipsoideGate contructors
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ellipsoideGate <- function(filterId="ellipsoidGate",.gate, distance,...) {
    if(missing(.gate) || !is.matrix(.gate))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
      
    new("ellipsoidGate",filterId=filterId,parameters=colnames(.gate), focus=.gate,distance=distance)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
