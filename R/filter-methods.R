## last changed November 2, 2006
## pdh: rectangleGate was still using cytoFrame. I changed it to flowFrame
##      made the filterDetails more uniform among the methods
##      got the right filename for a flowFrame to go in the msg

## ==========================================================================

## This implementation of filter does not require a per-filter implementation
## anything filter specific has been moved to either filterDetails or the %in%
## implementation. A filter operation is now about assembling those two calls 
## AND associating a filterResult with a particular flowFrame if possible.
setMethod("filter",signature("flowFrame","filter"),function(x,filter) {
	result                = as(x %in% filter,"filterResult")
	filterDetails(result) = filter
	result@frameId        = identifier(x)
	result
})
setReplaceMethod("filterDetails",signature("filterResult","filter"),function(result,value) {
	result@parameters    = value@parameters
	result@filterId      = value@filterId
	result@filterDetails = list(filter=value)
	
	result
})

#msg = paste(class(filter),"applied on",
#	deparse(substitute(flowObject))," (file:",
#	basename(flowObject@description["$FIL"]),") a",
#	class(flowObject),"object")


## Printing out a filter should give us something at least mildly sensible.
setMethod("show","filter",function(object) 
	cat(paste("A filter named '",object@filterId,"'\n",sep="")))
## Most filters define only a single population
setMethod("length","filter",function(x) 1)


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
        ##If there is only a single dimension then we have a degenerate case.
	if(ndim==1) 
		!is.na(cut(exprs(x)[,table@parameters[[1]]],range(table@boundaries[,1]),labels=FALSE))
	else if(ndim==2) {
            print(table@boundaries)
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
			if(nrow(y)>table@n) covMcd(y[sample(nrow(y),table@n),]) else covMcd(y)
		},
		cov.rob={cov.rob(y)},
		stop("How did you get here?")
		)
	W  = t(y)-cov$center
	result = exp(-.5*colSums((solve(cov$cov)%*%W)*W))>exp(-.5*table@scale.factor^2)
	attr(result,'center') = cov$center
	attr(result,'cov')    = cov$cov
	result
})
setReplaceMethod("filterDetails",signature("filterResult","norm2Filter"),function(result,value) {
	result = callNextMethod()
	#Bring the cov and center values into the details for easier access
	result@filterDetails$cov = attr(result@subSet,'cov')
	result@filterDetails$center = attr(result@subSet,'center')
	result
})

setMethod("summary","filter",function(object,...) {
	l = as(object,"logical")
	true = sum(l)
	count= length(l)
	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),name=object@filterId),class="filterSummary")
})
setMethod("summary","subsetFilter",function(object,...) {
	e1 = as(object@left,"logical")
	e2 = as(object@right,"logical")
	true = sum(e1&e2)
	count= sum(e2)
	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),name=object@filterId),class="filterSummary")	
})
