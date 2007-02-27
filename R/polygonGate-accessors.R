## ==========================================================================
## R wrapper for C function inPolygon
## checks for input arguments and makes sure that the polygon is closed
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inpolygon <- function(points, vertices){
  ## check validity of arguments
  dp <- dim(points)
  if(!is.matrix(points) || dp[1]<1 | dp[2]!=2)
    stop("Argument 'points' must be numeric matrix of two columns and at least one row ",
         "specifiying points on a two-dimensional plane")
  dv <- dim(vertices)
  if(!is.matrix(vertices) || dv[1]<2 | dv[2]!=2)
    stop("Argument 'vertices' must be numeric matrix of two columns and at least two rows",
         " specifying vertices of a polygon on a two-dimensional plane")
  
  ## the polygon must be closed
  if(!all(vertices[1,] == vertices[dv[1],]))
    vertices <- rbind(vertices, vertices[1,])

  ## call C function
  .Call("inpolygon", points, vertices, package="flowCore")
}




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
               
	 	as.logical(flowCore:::inPolygon(exprs(x)[,table@parameters],
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
