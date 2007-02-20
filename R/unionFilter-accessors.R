## ==========================================================================
## A union filter returns TRUE if ANY of the argument filters return true.
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","unionFilter"),function(x,table) {	
	fr = sapply(table@filters,filter,x=x)
	res= apply(sapply(fr,as,"logical"),1,any)
	details = list()
	for(i in fr) { 
		fd = filterDetails(i)
		details[names(fd)] = fd
	}
	attr(res,'filterDetails') = details
	res
})
