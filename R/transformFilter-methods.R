## Constructors
setMethod("%on%",signature("filter","transformList"),function(e1,e2) {
	new("transformFilter",
		filterId=paste(e1@filterId,"on transformed values of",paste(colnames(e2),collapse=",")),
		transforms=e2,
		filter=e1,parameters=e1@parameters)
})

# Not needed
#setMethod("summarizeFilter",signature("filterResult","transformFilter"),function(result,filter) {
#	ret = callNextMethod()
#	
#})

setMethod("%in%",signature("flowFrame","transformFilter"),function(x,table) {
	(table@transforms %on% x) %in% table@filter
})
