## Constructors
kmeansFilter = function(...,filterId="kmeans") {
	l = length(list(...))
	if(l>1)
		stop("Barcode filters only operate on a single parameter.")
	x = ..1
	if(is.list(x)) {
	} else
		new("kmeansFilter",parameters=names(list(...))[1],populations=x,filterId=filterId)
}

## Filtering Methods -- we are not a logical filter so we return a vector
## of indices indicating a population.
setMethod("filterDetails",signature("flowFrame","kmeansFilter","ANY"),function(flowObject,filter,result) {
	msg = paste(class(filter),"applied on",
		deparse(substitute(flowObject))," (file:",
		basename(flowObject@description["$FIL"]),") a",
		class(flowObject),"object")
	list(message=msg,filter=filter,populations=filter@populations)	
})
setMethod("%in%",signature("flowFrame","kmeansFilter"),function(x,table) {
	##We accomplish the actual filtering via K-means
	param  = table@parameters[1]
	values = exprs(x)[,param]
	npop   = length(table@populations)
	km     = kmeans(values,centers=quantile(values,(1:npop)/(npop+1)))
	#Ensure that the populations are sorted according to their center,
	#which matches the assumption of the population vector.
	structure(as.integer(order(km$centers)[km$cluster]),levels=table@populations)
})
setMethod("length",signature("kmeansFilter"),function(x) length(x@populations))
setMethod("show",signature("kmeansFilter"),function(object) {
	msg = paste("A k-means filter named '",object@filterId,"': ",table@parameters[1]," (",
		paste(table@populations,sep=","),")")
	cat(msg)
	cat("\n")
	invisible(msg)
})
