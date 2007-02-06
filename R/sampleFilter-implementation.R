## Constructors
sampleFilter = function(filterId="sample",size) {
	new("sampleFilter",parameters=character(0),filterId=filterId,size=size)
}
setMethod("%in%",signature("flowFrame","sampleFilter"),function(x,table) {
	n = if(table@size > nrow(x)) nrow(x) else table@size
	l = rep(FALSE,nrow(x))
	l[sample(length(l),n,replace=FALSE)] = TRUE
	l
})
setMethod("show",signature("sampleFilter"),function(object) {
	msg = paste("A filter named '",object@filterId,"' returning objects with ",object@size," rows",sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})
