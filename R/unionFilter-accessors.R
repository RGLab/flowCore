## ==========================================================================
## A union filter returns TRUE if ANY of the argument filters return true.
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","unionFilter"),function(x,table) {	
	apply(sapply(table@filters,"%in%",x=x),1,any)
})
