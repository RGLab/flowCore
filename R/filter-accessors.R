## ==========================================================================
## This implementation of filter does not require a per-filter implementation
## anything filter specific has been moved to either filterDetails or the %in%
## implementation. A filter operation is now about assembling those two calls 
## AND associating a filterResult with a particular flowFrame if possible.
## ==========================================================================


## ==========================================================================
## accessor method for slot parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",signature("filter"),function(object) character(0))
setMethod("parameters",signature("parameterFilter"),function(object) object@parameters)
setMethod("parameters",signature("setOperationFilter"),function(object) unique(unlist(lapply(object@filters,parameters))))

## ==========================================================================
## summary method for filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",signature("filter"),function(object,result,...) {
	l = if(missing(result)) as(object,"logical") else as(result,"logical")
	true=sum(l)
	count=length(l)
	new("filterSummary",name=identifier(object),true=true,count=count,p=c/l)
})

# setMethod("summary",signature("filter"),function(object,result,...) {
# 	l = if(missing(result)) as(object,"logical") else as(result,"logical")
# 	true = sum(l)
# 	count= length(l)
# 	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),
#                        name=identifier(object)),class="filterSummary")
# })


## Printing out a filter should give us something at least mildly sensible.
setMethod("show",signature("filter"),function(object) 
	cat(paste("A filter named '",object@filterId,"'\n",sep="")))

## Most filters define only a single population
## setMethod("length","filter",function(x) 1)
setMethod("identifier", signature="filter",
          definition=function(object) object@filterId)
setReplaceMethod("identifier",signature("filter","character"),function(object,value) {
	object@filterId = value
	object
})


## ==========================================================================
## The %subset% operator constructs a filter that takes the subset of the LHS
## filter in the RHS filter result. 
## For many cases this is equivalent to an intersection filter.
## --------------------------------------------------------------------------
setMethod("%subset%",signature("filter","filter"),
          function(e1,e2) {
            new("subsetFilter",
                filters=list(e1,e2),filterId=paste(identifier(e1),"in",identifier(e2)))
})
setMethod("%subset%",signature("list","filter"),function(e1,e2) lapply(e1,"%subset%",e2=e2))
setMethod("%subset%",signature("filterSet","filter"),function(e1,e2) {
	#Make a copy of the filterSet, preserving R semantics
	x = as(as(e1,"list"),"filterSet")
	n = names(e1)
	x[[""]] = e2
	target = as.symbol(identifier(e2))
	for(i in n) {
		x[[""]] = substitute(~ a %subset% b,list(a=as.symbol(i),b=target))
	}
	x
})

## ==========================================================================
## The | operator returns an object representing the union of two filters.
## A simple optimization would be to linearize the union of a filter and
## another union filter.
## --------------------------------------------------------------------------
setMethod("|",signature("filter","filter"),
          function(e1,e2) {
            new("unionFilter",
                filters=list(e1,e2),filterId=paste(identifier(e1),"or",identifier(e2)))
          })
setMethod("|",signature("list","filter"),function(e1,e2) lapply(e1,"|",e2=e2))
setMethod("|",signature("filter","list"),function(e1,e2) lapply(e2,"|",e1=e1))

## ==========================================================================
## The & operator returns an object representing the intersection of two
## filters. This is somewhat different from the %subset% operator because
## some filters depend on the  data and would return different results
## when applied to the full dataset.
## --------------------------------------------------------------------------
setMethod("&",signature("filter","filter"),
          function(e1,e2) {
            new("intersectFilter",
                filters=list(e1,e2),filterId=paste(identifier(e1),"and",identifier(e2)))
          })
setMethod("&",signature("list","filter"),function(e1,e2) lapply(e1,"&",e2=e2))
setMethod("&",signature("filter","list"),function(e1,e2) lapply(e2,"&",e1=e1))



## ==========================================================================
## The ! operator returns an object that simply takes the complement of the
## filter it returns.
## --------------------------------------------------------------------------
setMethod("!",signature("filter"),function(x) {
	new("complementFilter",filters=list(x),
            filterId=paste("not",identifier(x)))
})


## ==========================================================================
## --------------------------------------------------------------------------
"%&%" = function(e1,e2) e1 %subset% e2

## ==========================================================================
## Plotting method for filters lets us get a basic plot of the data and a 
## filter. Uses the draw method to actually display the filter on some data.
## --------------------------------------------------------------------------
setMethod("plot",signature(x="flowFrame", y="filter"),function(x, y, z=NULL, results=NULL, ...) {
    if(is.null(z))	z = parameters(y)
    if(length(z)>2) {
        warning("Only plotting the first two parameters.")
        z = z[1:2]
    }
    ##First we plot the actual data.
    plot(x, z, ...)
    ##Then draw the filter on top of it which takes the form (filter, data, params, results (if any))
    draw(y, x, z, results, ...)
})

setMethod("draw", signature("filter","flowFrame"), function(x, data, params=NULL, results=NULL, ...) {
	if(is.null(params))
          params = parameters(x)
	if(length(params)>2) {
		params = params[1:2]
	}
        if(is(x,"rectangleGate")){
          lines(c(x@min[params[2]], x@min[params[2]]), c(x@min[params[1]], x@max[params[1]]), col="red", ...)
          lines(c(x@min[params[1]], x@min[params[1]]), c(x@min[params[2]], x@max[params[2]]), col="red", ...)
      }
        if(is(x,"polygonGate") | is(x,"polytopeGate"))
          lines(x@boundaries[,params[1]], x@boundaries[,params[2]], col="red", ...)
        ##if(is(x, ellipsoidGate))
          
})
