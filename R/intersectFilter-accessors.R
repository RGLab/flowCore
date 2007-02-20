## ==========================================================================
## An intersectFilter only returns TRUE if ALL the member filters are TRUE.
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","intersectFilter"),function(x,table) {
	fr = sapply(table@filters,filter,x=x)
	res= apply(sapply(fr,as,"logical"),1,all)
	details = list()
	for(i in fr) { 
		fd = filterDetails(i)
		details[names(fd)] = fd
	}
	attr(res,'filterDetails') = details
	res})
