
## ==========================================================================
## The | operator returns an object representing the union of two filters. A simple
## optimization would be to linearize the union of a filter and another union filter.
setMethod("|",signature("filter","filter"),function(e1,e2) {
	new("unionFilter",parameters=unique(c(e1@parameters,e2@parameters)),filters=list(e1,e2),filterId=paste(e1@filterId,"or",e2@filterId))
})
## A union filter returns TRUE if ANY of the argument filters return true.
setMethod("%in%",c("flowFrame","unionFilter"),function(x,table) {	
	apply(sapply(table@filters,"%in%",x=x),1,any)
})

## The & operator returns an object representing the intersection of two filters. This
## is somewhat different from the %subset% operator because some filters depend on the 
## data and would return different results when applied to the full dataset.
setMethod("&",signature("filter","filter"),function(e1,e2) {
	new("intersectFilter",parameters=unique(c(e1@parameters,e2@parameters)),filters=list(e1,e2),filterId=paste(e1@filterId,"and",e2@filterId))
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
	new("subsetFilter",parameters=unique(c(e1@parameters,e2@parameters)),left=e1,right=e2,filterId=paste(e1@filterId,"in",e2@filterId))
})
setMethod("%&%",signature("ANY","ANY"),function(e1,e2) e1 %subset% e2)
## Returns TRUE for elements that are true on the LHS and the RHS, however the LHS filter is only executed against the subset returned
## by the RHS filtering operation. This is particularly important for unsupervised filters like norm2Filter. The result is still relative
## to the ENTIRE flowFrame however. 
setMethod("%in%",c("flowFrame","subsetFilter"),function(x,table) {
	y = filter(x,table@right)
	if(length(y) == 1) {
		y = as(y,"logical")
		n = which(y)
		y[n[!(x[n,] %in% table@left)]] = FALSE
		y
	} else if(length(y) > 1){
		res = rep(NA,nrow(x))
		for(i in seq(along=y)) {
			z = as(y[[i]],"logical")
			n = which(z)
			z[n[!(x[n,] %in% table@left)]] = FALSE
			res[z] = i
		}
		structure(as.integer(res),levels=paste(identifier(table@left),"in",names(y)))
	}
})
setMethod("%in%",c("flowFrame","filterResult"),function(x,table) {
	frameId = identifier(x)
	if(all(!is.na(c(frameId,table@frameId))) && frameId != table@frameId)
		warning("Frame identifiers do not match. It is possible that this filter is not compatible with this frame.")
	if(nrow(x) != length(table@subSet))
		stop("Number of rows in frame do not match those expected by this filter.")
	table@subSet
})

#Allow the coercion of resolvable filters (i.e. those derived from filterResult) to be composed and then
#converted into a logical vector. This allows for a lot of processing to be done simply using the filter
#results.
setAs("filter","logical",function(from) stop("Only resolved filters can be converted to a logical vector."))
setAs("subsetFilter","logical",function(from) as(from@left,"logical") & as(from@right,"logical"))
setAs("intersectFilter","logical",function(from) apply(sapply(from@filters,as,Class="logical"),1,all))
setAs("unionFilter","logical",function(from) apply(sapply(from@filters,as,Class="logical"),1,any))
setAs("complementFilter","logical",function(from) !as(from@filter,"logical"))
