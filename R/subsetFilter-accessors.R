## ==========================================================================
## Methods for objects of type 'intersectFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("subsetFilter"),
          function(object)
      {
          cat("filter '", identifier(object),
              "'\nthe filtering operation defined by\n", sep="")
          print(object@filters[[1]])
          cat("after subsetting by\n")
          print(object@filters[[2]])
      })



# FIXME: What do summarizeFilter methods do and why are they not
## exported in the name space?
setMethod("summarizeFilter",signature("filterResult","subsetFilter"),function(result,filter) {
	ret = callNextMethod()
	ret$subsetCount = attr(result@subSet,'subsetCount')
	ret
})


## ==========================================================================
## summary method for subsetFilter
## --------------------------------------------------------------------------
setMethod("summary","subsetFilter",function(object,result,...) {
	if(missing(result)) {
		e1 = as(object@filters[[1]],"logical")
		e2 = as(object@filters[[2]],"logical")
		true = sum(e1&e2)
		count= sum(e2)		
		new("filterSummary",name=identifier(object),true=true,count=count,p=true/count)			
	} else {
		true = sum(as(result,"logical"))
		count = filterDetails(result,identifier(object))$subsetCount
		new("filterSummary",name=identifier(object),true=true,count=count,p=true/count)			
	}
})
