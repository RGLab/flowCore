## ==========================================================================
## %in% method
## ---------------------------------------------------------------------------
setMethod("%in%",signature("flowFrame","sampleFilter"),function(x,table) {
	n = if(table@size > nrow(x)) nrow(x) else table@size
	l = rep(FALSE,nrow(x))
	l[sample(length(l),n,replace=FALSE)] = TRUE
	l
})


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("sampleFilter"),function(object) {
	msg = paste("A filter named '",object@filterId,"' returning objects with ",object@size," rows",sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})
