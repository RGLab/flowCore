## Constructors
kmeansFilter = function(...,filterId="kmeans") {
	l = length(list(...))
	if(l>1)
		stop("k-means filters only operate on a single parameter.")
	x = ..1
	if(is.list(x)) {
		new("kmeansFilter",parameters=names(x)[1],populations=x[[1]],filterId=filterId)
	} else
		new("kmeansFilter",parameters=names(list(...))[1],populations=x,filterId=filterId)
}

## Filtering Methods -- we are not a logical filter so we return a vector
## of indices indicating a population.
setMethod("summarizeFilter",signature("filterResult","kmeansFilter"),function(result,filter) {
	ret = callNextMethod()
	ret$populations = filter@populations
	ret
})
setMethod("%in%",signature("flowFrame","kmeansFilter"),function(x,table) {
	##We accomplish the actual filtering via K-means
	param  = table@parameters[1]
	values = exprs(x)[,param]
	npop   = length(table@populations)
	km     = kmeans(values,centers=quantile(values,(1:npop)/(npop+1)))
	#Ensure that the populations are sorted according to their center,
	#which matches the assumption of the population vector.
	structure(as.integer(order(km$centers)[km$cluster]),class="factor",levels=table@populations)
})
setMethod("length",signature("kmeansFilter"),function(x) length(x@populations))
setMethod("show",signature("kmeansFilter"),function(object) {
	msg = paste("A k-means filter named '",object@filterId,"': ",object@parameters[1]," (",
		paste(object@populations,collapse=","),")")
	cat(msg)
	cat("\n")
	invisible(msg)
})
