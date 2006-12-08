## ==========================================================================
## Define population inclusion methods 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setAs("filter","population",function(from) new("filterPopulation",filter=from))

# populations act a lot like filters
setMethod("%in%",signature(x="ANY",table="filterPopulation"),function(x,table) x %in% table@filter)
setMethod("%in%",signature(x="ANY",table="unionPopulation"),function(x,table) apply(sapply(table@members,"%in%",x=x),1,any))
setMethod("%in%",signature(x="ANY",table="intersectionPopulation"),function(x,table) apply(sapply(table@members,"%in%",x=x),1,all))
setMethod("%in%",signature(x="ANY",table="complementPopulation"),function(x,table) !(x %in% table@population))
## ==========================================================================



## ==========================================================================
## Compose populations using logic filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("|",signature("population","population"),function(e1,e2) new("unionPopulation",members=list(e1,e2)))
setMethod("&",signature("population","population"),function(e1,e2) new("intersectionPopulation",members=list(e1,e2)))
setMethod("!",signature("population"),function(e1) new("complementPopulation",population=e1))
## ==========================================================================


## ==========================================================================
## show methods for for Population
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show","filterPopulation",function(object) {
	cat("A population defined by the filter:\n")
	show(object@filter)
})
setMethod("show","complementPopulation",function(object) {
	cat("A population defined by the complement of the following population:\n")
	show(object@population)
})
setMethod("show","unionPopulation",function(object) {
	cat("A population defined by the union of the following populations:\n")
	for(i in object@members) show(i)
})
setMethod("show","intersectionPopulation",function(object) {
	cat("A population defined by the intersection of the following populations:\n")
	for(i in object@members) show(i)
})
## ==========================================================================
