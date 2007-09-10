#Everything is just passed to the referenced filter. This may be better handled by the type system by
#having "real" filters inherit from concreteFilter (or something) and then simply having a setAs(), but 
#I think that will be too much work for filter authors.
setMethod("identifier",signature("filterReference"),function(object) {
	if(exists(object@name,env=object@env)) identifier(as(object,"concreteFilter")) else object@name
})

setMethod("parameters",signature("filterReference"),
	function(object) parameters(as(object,"concreteFilter")))
setMethod("summary",signature("filterReference"),
	function(object,...) summary(as(object,"concreteFilter"),...))
setMethod("show",signature("filterReference"),function(object) {
	if(exists(object@name,env=object@env))
		cat(paste("A reference to a filter named '",identifier(object),"'\n",sep=""))
	else
		cat(paste("An unresolvable reference to a filter named '",object@name,"'\n",sep=""))
})
setMethod("%in%",signature("ANY","filterReference"),function(x,table) x %in% as(table,"concreteFilter"))
setMethod("summarizeFilter",signature("filterResult","filterReference"),
	function(result,filter) summarizeFilter(result,as(filter,"concreteFilter")))

setMethod("filterReference",signature("environment","character"),function(from,name) {
	new("filterReference",name=name,env=from)
})