## last changed November 2, 2006
## pdh: rectangleGate was still using cytoFrame. I changed it to flowFrame
##      made the filterDetails more uniform among the methods
##      got the right filename for a flowFrame to go in the msg

## ==========================================================================

## This implementation of filter does not require a per-filter implementation
## anything filter specific has been moved to either filterDetails or the %in%
## implementation. A filter operation is now about assembling those two calls 
## AND associating a filterResult with a particular flowFrame if possible.
setMethod("filter",signature("flowFrame","filter"),function(flowObject,filter) {
	result = flowObject %in% filter
	details= filterDetails(flowObject,result,filter)
	new("filterResult",
		parameters=filter@parameters,
		filterId=filter@filterId,
		frameId=if(is.null(flowObject@description["GUID"])) flowObject@description["$FIL"] else flowObject@description["GUID"],
		subSet=result,filterDetails=detauls)
})
## Printing out a filter should give us something at least mildly sensible.
setMethod("show","filter",function(object) 
	cat(paste("A filter named '",object@filterId,"'\n",sep="")))
## ==========================================================================
## filterDetails constructs the metadata carried around by the filterResult
## this should generally be specific enough, but you may want to override your
## own.
setMethod("filterDetails",signature("flowFrame","filter","ANY"),function(flowObject,filter,result) {
	msg = paste(class(filter),"applied on",
		deparse(substitute(flowObject))," (file:",
		basename(flowObject@description["$FIL"]),") a",
		class(flowObject),"object")
	list(message=msg,filter=filter)
})


## ==========================================================================
## support methods for filtering flowFrame Object using rectangleGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature(x="flowFrame",table="rectangleGate"),function(x,table) {
	e = if(length(table@parameters)==1) as.matrix(exprs(x)[,table@parameters]) else exprs(x)[,table@parameters]
	apply(sapply(seq(along=table@parameters),function(i) {
		!is.na(cut(e[,i],c(table@min[i],table@max[i]),labels=FALSE))
	}),1,all)
})
setMethod("show","rectangleGate",function(object) {
	cat("Rectangular gate with dimensions:\n")
	for(i in seq(along=object@parameters)) {
		cat("\t")
		cat(object@parameters[i])
		cat(": (")
		cat(paste(object@min[i],object@max[i],sep=","))
		cat(")\n")
	}
})

#setMethod("filter",
#          signature=signature(filter="rectangleGate",flowObject="flowFrame",parent="ANY"),
#          definition=function(filter,flowObject,parent) {
#              selectNew <-rectangleFiltering(filter,flowObject,parent)
#               msg <- paste("rectangleGate applied on ",
#                            deparse(substitute(flowObject)),
#                            " (file:",basename(flowObject@description["$FIL"]),
#                            ") a ", class(flowObject)," object", sep="")
#               out = list(msg,filter=filter)
#              new("filterResult", subSet=selectNew,filterDetails=out)
#          })
## ==========================================================================



## ==========================================================================
## filter flowFrame Object using polygonGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%in%",signature("flowFrame","polygonGate"),function(x,table) {
	ndim = length(table@parameters)
	#If there is only a single dimension then we have a degenerate case.
	if(ndim==1) 
		!is.na(cut(exprs(x)[,table@parameters[[1]]],range(table@boundaries[,1]),labels=FALSE))
	else if(ndim==2) {
	 	as.logical(.Call(inPolygon,exprs(x)[,table@parameters],table@boundaries))
	} else 
		stop("Polygonal gates only support 1 or 2 dimensional gates (for now).")
})


#setMethod("filter",
#          signature=signature(filter="polygonGate",flowObject="flowFrame",parent="ANY"),
#          definition=function(flowObject,filter,parent) {
#              selectNew <-polygonFiltering(filter,flowObject,parent)
#              msg <- paste("polygonGate applied on ",
#                           deparse(substitute(flowObject)),
#                           " (file:",basename(flowObject@description["$FIL"]),
#                           ") a ", class(flowObject)," object", sep="")
#       			out = list(msg,filter=filter)
#              new("filterResult", subSet=selectNew,filterDetails=out)
#          } 
#          )


## ==========================================================================

## ==========================================================================
##  filter  flowFrame Object using norm2Filter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#setMethod("filter",
#          signature=signature(filter="norm2Filter",flowObject="flowFrame",parent="ANY"),
#          definition=function(filter,flowObject,parent) {
#              selectNew <-normFiltering(filter,flowObject,parent)
#              msg <- paste("norm2Filter applied on ",
#                            deparse(substitute(flowObject)),
#                            " (file:",basename(flowObject@description["$FIL"]),
#                            ") a ", class(flowObject)," object", sep="")
#               out = list(msg,filter=filter,mu=selectNew$mu,S=selectNew$S)
#              new("filterResult", subSet=selectNew$sel,filterDetails=out)
#          })
## ==========================================================================
setMethod("%in%",signature("flowFrame",table="norm2Filter"),function(x,table) {
	if(length(table@parameters) != 2)
		stop("norm2 filters require exactly two parameters.")
	y = if(length(table@transformation)>0) 
		exprs(do.call("transform",c(x,table@transformation)))[,table@parameters] else exprs(x)[,table@parameters]
	
	if(is.na(match(table@method,c("covMcd","cov.rob"))))
		stop("Method must be either 'covMcd' or 'cov.rob'")
	cov = switch(table@method,
		covMcd = {
			if(nrow(y)>50000) covMcd(y[sample(50000,nrow(y)),]) else covMcd(y)
		},
		cov.rob={cov.rob(y)},
		stop("How did you get here?")
		)
	W  = t(y)-cov$center
	exp(-.5*colSums((solve(cov$cov)%*%W)*W))>exp(-.5*table@scale.factor^2)
})

## ==========================================================================
## The | operator returns an object representing the union of two filters. A simple
## optimization would be to linearize the union of a filter and another union filter.
setMethod("|",signature("filter","filter"),function(e1,e2) {
	p1 = sapply(e1@parameters,"grep",x=e2@parameters)	
	new("unionFilter",parameters=c(e1@parameters[-p1],e2@parameters),filters=list(e1,e2),filterId=paste(e1@filterId,"or",e2@filterId))
})
## A union filter returns TRUE if ANY of the argument filters return true.
setMethod("%in%",c("flowFrame","unionFilter"),function(x,table) {	
	apply(sapply(table@filters,"%in%",x=x),1,any)
})

## The & operator returns an object representing the intersection of two filters. This
## is somewhat different from the %subset% operator because some filters depend on the 
## data and would return different results when applied to the full dataset.
setMethod("&",signature("filter","filter"),function(e1,e2) {
	p1 = sapply(e1@parameters,"grep",x=e2@parameters)
	new("intersectFilter",parameters=c(e1@parameters[-p1],e2@parameters),filters=list(e1,e2),filterId=paste(e1@filterId,"and",e2@filterId))
})
## An intersectFilter only returns TRUE if ALL the member filters are TRUE.
setMethod("%in%",c("flowFrame","intersectFilter"),function(x,table) {
	apply(sapply(table@filters,"%in%",x=x),1,all)
})

## The ! operator returns an object that simply takes the complement of the filter it returns.
setMethod("!",signature("filter"),function(e1) 
	new("complementFilter",parameters=e1@parameters,filter=e1,filterId=paste("not",e1@filterId)))
## Returns TRUE when the input filter is FALSE.
setMethod("%in%",c("flowFrame","complementFilter"),function(x,table) !(x%in%table@filter))

## The %subset% operator constructs a filter that takes the subset of the LHS filter in the RHS filter result. 
## For many cases this is equivalent to an intersection filter.
setMethod("%subset%",signature("filter","filter"),function(e1,e2) {
	p1 = sapply(e1@parameters,"grep",x=e2@parameters)
	new("subsetFilter",parameters=c(e1@parameters[-p1],e2@parameters),left=e1,right=e2,filterId=paste(e1@filterId,"in",e2@filterId))
})
## Returns TRUE for elements that are true on the LHS and the RHS, however the LHS filter is only executed against the subset returned
## by the RHS filtering operation. This is particularly important for unsupervised filters like norm2Filter. The result is still relative
## to the ENTIRE flowFrame however. 
setMethod("%in%",c("flowFrame","subsetFilter"),function(x,table) {
	y = x %in% table@right
	n = which(y)
	y[n[!(x[n,] %in% table@left)]] = FALSE
	y
})

setMethod("%in%",c("flowFrame","filterResult"),function(x,table) {
	frameId = if(is.null(x@description["GUID"])) x@description["$FIL"] else x@description["GUID"]
	if(frameId != table@frameId)
		warning("Frame identifiers do not match. It is possible that this filter is not compatible with this frame.")
	if(nrow(x) != length(table@subSet))
		stop("Number of rows in frame do not match those expected by this filter.")
	table@subSet
})

#Allow the coercion of resolvable filters (i.e. those derived from filterResult) to be composed and then
#converted into a logical vector. This allows for a lot of processing to be done simply using the filter
#results.
setAs("filter","logical",function(from) stop("Only resolved filters can be converted to a logical vector."))
setAs("filterResult","logical",function(from) if(is.logical(from@subSet)) from@subSet else from@subSet==1)
setAs("subsetFilter","logical",function(from) as(from@left,"logical") & as(from@right,"logical"))
setAs("intersectFilter","logical",function(from) apply(sapply(table@filters,as,Class="logical"),1,all))
setAs("unionFilter","logical",function(from) apply(sapply(table@filters,as,Class="logical"),1,any))
setAs("complementFilter","logical",function(from) !as(from@filter,"logical"))
