




## ==========================================================================
## length method
## ---------------------------------------------------------------------------
setMethod("length","manyFilterResult",function(x) ncol(x@subSet))


## ==========================================================================
## names method
## ---------------------------------------------------------------------------
setMethod("names","manyFilterResult",function(x) colnames(x@subSet))


## ==========================================================================
## subsetting method
## ---------------------------------------------------------------------------
setMethod("[[",signature("manyFilterResult"),function(x,i,j,drop=FALSE) {
	if(is.numeric(i)) i = names(x)[i]
	if(length(i)!=1)
		stop("Only a single subpopulation can be selected.")
	filterDetails = structure(list(filterDetails(x,i)),names=i)
	filterDetails$population = i
	filterDetails$source     = identifier(x)
	new("logicalFilterResult",subSet=x@subSet[,i],filterDetails=filterDetails,
		frameId=x@frameId,filterId=i)
})


as.data.frame.manyFilterResult = function(x,row.names=NULL,optional=FALSE,...) {
	nrows = length(x)
	nm    = if(nzchar(identifier(x))) identifier(x) else "filter"
	if(is.null(row.names)) {
		if(nrows == 0L)
			row.names = character(0L)
		else if(length(row.names <- names(x)) == nrows && !any(duplicated(row.names))) {
		} else row.names = .set_row_names(nrows)
	}
	values = sapply(names(x),function(i) {
		m = x[[i]]
		paste(deparse(as(filterDetails(m,identifier(m))$filter,"call")),sep="\n",collapse="\n")
	})
	names(values) = NULL
	values = list(values)
	if(!optional)
		names(values) = nm
	attr(values,"row.names") = row.names
	class(values) = "data.frame"
	values
}
