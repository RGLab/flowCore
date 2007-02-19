## ==========================================================================
## filter flowFrame object using polygonGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature(x="flowFrame",table="polygonGate"),function(x,table) {
	ndim = length(table@parameters)
        ##If there is only a single dimension then we have a degenerate case.
	if(ndim==1) 
		!is.na(cut(exprs(x)[,table@parameters[[1]]],
                           range(table@boundaries[,1]),labels=FALSE, right=FALSE))
	else if(ndim==2) {
               
	 	as.logical(.Call(inPolygon,exprs(x)[,table@parameters],
                                 table@boundaries))
            
	} else 
		stop("Polygonal gates only support 1 or 2 dimensional gates (for now).")
})


## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature(object="polygonGate"),function(object) {
  nb <-  nrow(object@boundaries)
  cat("Polygonal gate with", ifelse(all(is.na(object@boundaries)), 0, nb),
      "vertices\n")
})
