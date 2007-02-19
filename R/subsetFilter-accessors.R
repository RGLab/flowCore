## ==========================================================================
## Returns TRUE for elements that are true on the LHS and the RHS, however the
## LHS filter is only executed against the subset returned by the RHS filtering
## operation. This is particularly important for unsupervised filters like
## norm2Filter. The result is still relative to the ENTIRE flowFrame however.
## ---------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","subsetFilter"),function(x,table) {
	y = filter(x,table@filters[[2]])
	if(length(y) == 1) {
		y = as(y,"logical")
		n = which(y)
		y[n[!(x[n,] %in% table@filters[[1]])]] = FALSE
		y
	} else if(length(y) > 1){
		res = rep(NA,nrow(x))
		for(i in seq(along=y)) {
			z = as(y[[i]],"logical")
			n = which(z)
			z[n[!(x[n,] %in% table@filters[[1]])]] = FALSE
			res[z] = i
		}
		structure(as.integer(res),levels=paste(identifier(table@filters[[1]]),
                                            "in",names(y)))
	}
})


## ==========================================================================
## summary method for subsetFiler
## --------------------------------------------------------------------------
setMethod("summary","subsetFilter",function(object,...) {
	e1 = as(object@filters[[1]],"logical")
	e2 = as(object@filters[[2]],"logical")
	true = sum(e1&e2)
	count= sum(e2)
	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),
                       name=object@filterId),class="filterSummary")	
})
