## ==========================================================================
## Returns TRUE when the input filter is FALSE.
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","complementFilter"),function(x,table) {
	r = filter(x,table@filters[[1]])
	z = !as(r,"logical")
	attr(z,'filterDetails') = filterDetails(r)
	z
})

