## ==========================================================================
## Add or Replacement method to complete or replace the results
## coming out of a filtering operation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("filterDetails",signature("filterResult","missing"),function(result,filterId) {
	result@filterDetails
})
setMethod("filterDetails",signature("filterResult","ANY"),function(result,filterId) {
	result@filterDetails[[filterId]]
})

setReplaceMethod("filterDetails",signature("filterResult","character",value="ANY"),
                 function(result,filterId,...,value) {
	result@filterDetails[[filterId]] = value
	result
})
setReplaceMethod("filterDetails",signature("filterResult","character",value="filter"),
                 function(result,filterId,...,value) {
	filterDetails(result,filterId) = summarizeFilter(result,value)
	result
})
setReplaceMethod("filterDetails",signature("filterResult","character",value="setOperationFilter"),
                 function(result,filterId,...,value) {
	details = attr(result@subSet,'filterDetails')
	for(i in names(details)) {
		filterDetails(result,i) = details[[i]]
	}
	#Record ourselves for posterity
	filterDetails(result,filterId) = summarizeFilter(result,value)
	result
})


## ==========================================================================
## Summarize a filtering operation
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summarizeFilter",signature("filterResult","filter"),function(result,filter) {
	list(filter=filter)
})
setMethod("summarizeFilter",signature("filterResult","parameterFilter"),function(result,filter) {
	ret = callNextMethod()
	ret$parameters = parameters(filter)
	ret
})

setMethod("summary",signature("filterResult"),function(object,...)
	summary(filterDetails(object,object@filterId)$filter,object,...))


## ==========================================================================
## --------------------------------------------------------------------------
setMethod("%in%",c("flowFrame","filterResult"),function(x,table) {
	as(table,"logical")
})

## ==========================================================================
## Allow us to compare filterResults and flowFrames. This lets us check
## (and warn or stop) that a particular flowFrame generated a filterResult
## allowing us to use it for further processing.
## --------------------------------------------------------------------------
setMethod("==",signature("flowFrame","filterResult"),
          definition=function(e1,e2) {
            i1 = identifer(e1)
            i2 = e2@frameId
            (length(i1) == 0 || length(i2) == 0 || i1 == i2)
          })
#Does S4 do this for us automagically? I don't know!
setMethod("==",signature("filterResult","flowFrame"),
          definition=function(e1,e2) e2==e1)

## ==========================================================================
## identifier method for filterResult
## --------------------------------------------------------------------------
setMethod("identifier", signature="filterResult",
          definition=function(object) object@filterId)


setMethod("[[",signature("filterResult"),function(x,i,j,drop=FALSE) {
	if((is.character(i) && i != x@filterId) || (as.numeric(i) > 1))
		stop("filter index out of bounds")
	x
})