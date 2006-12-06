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

rectGate = function(.gate,...,id="rectangleGate") {
	if(missing(.gate) || !is.matrix(.gate))
	.gate = sapply(if(missing(.gate)) list(...) else .gate,function(x) c("min"=x[1],"max"=x[2]))
	new("rectangleGate",filterId=id,parameters=colnames(.gate),min=.gate[1,],max=.gate[2,])
}



#modeGate = function(.gate,...,id="modeGate")=======

