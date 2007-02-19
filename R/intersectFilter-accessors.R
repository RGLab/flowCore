## ==========================================================================
## An intersectFilter only returns TRUE if ALL the member filters are TRUE.
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","intersectFilter"),function(x,table) {
	apply(sapply(table@filters,"%in%",x=x),1,all)
})
