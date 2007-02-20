## ==========================================================================
## Returns TRUE for elements that are true on the LHS and the RHS, however the
## LHS filter is only executed against the subset returned by the RHS filtering
## operation. This is particularly important for unsupervised filters like
## norm2Filter. The result is still relative to the ENTIRE flowFrame however.
## ---------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","subsetFilter"),function(x,table) {
	y = filter(x,table@filters[[2]])
	z = if(length(y) == 1) {
			w = as(y,"logical")
			n = which(w)
			r = filter(x[n,],table@filters[[1]])
			filterDetails(y,identifier(table@filters[[1]])) = summarizeFilter(r,table@filters[[1]])
			w[n[!as(r,"logical")]] = FALSE
			w
		} else {
			res = rep(NA,nrow(x))
			ll  = paste(identifier(table@filters[[1]]),"in",names(y))
			for(i in seq(along=y)) {
				w = as(y[[i]],"logical")
				n = which(w)
				r = filter(x[n,],table@filters[[1]])
				filterDetails(y,ll[i]) = summarizeFilter(r,table@filters[[1]])
				w[n[!as(r,"logical")]] = FALSE
				res[z] = i
			}
			structure(as.integer(res),levels=ll)
		}
	#We need to track our filterDetails to a higher level for summarizeFilter.
	attr(z,'filterDetails') = filterDetails(y)
	z
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
