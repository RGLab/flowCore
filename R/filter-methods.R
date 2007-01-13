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
	result@filterId       = filter@filterId
	result@parameters     = filter@parameters
	filterDetails(result,result@filterId) = filter
	result@frameId        = identifier(x)
	result
})
setReplaceMethod("filterDetails",signature("filterResult","character",value="ANY"),function(result,filterId,...,value) {
	result@filterDetails[[filterId]] = value
	result
})
setReplaceMethod("filterDetails",signature("filterResult","character",value="filter"),function(result,filterId,...,value) {
	filterDetails(result,filterId) = summarizeFilter(result,value)
	result
})

setMethod("filterDetails",signature("filterResult","missing"),function(result,filterId) {
	result@filterDetails
})
setMethod("filterDetails",signature("filterResult","ANY"),function(result,filterId) {
	result@filterDetails[[filterId]]
})


setMethod("summarizeFilter",signature("filterResult","filter"),function(result,filter) {
	list(filter=filter)
})

setMethod("summary","filter",function(object,...) {
	l = as(object,"logical")
	true = sum(l)
	count= length(l)
	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),name=object@filterId),class="filterSummary")
})

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
setMethod("summary","rectangleGate",function(object,...) {
	structure(sapply(seq(along=object@parameters),function(i) c(min=object@min[i],max=object@max[i])),names=object@parameters)
})



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

## ==========================================================================
##  filter  flowFrame Object using norm2Filter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
	attr(result,'radius') = table@scale.factor
	result
})
setMethod("summarizeFilter",signature("filterResult","norm2Filter"),function(result,filter) {
	ret = callNextMethod()
	ret$cov    = attr(result@subSet,'cov')
	ret$center = attr(result@subSet,'center')
	ret$radius = attr(result@subSet,'radius')
	ret
})

