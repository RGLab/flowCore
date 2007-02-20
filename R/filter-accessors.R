## ==========================================================================
## This implementation of filter does not require a per-filter implementation
## anything filter specific has been moved to either filterDetails or the %in%
## implementation. A filter operation is now about assembling those two calls 
## AND associating a filterResult with a particular flowFrame if possible.
## ==========================================================================


## ==========================================================================
## summary method for filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",signature("filter"),function(object,...) {
	l = as(object,"logical")
	true = sum(l)
	count= length(l)
	structure(list(true=true,false=count-true,n=count,p=true/count,q=1-(true/count),
                       name=object@filterId),class="filterSummary")
})


## ==========================================================================
## Summarize the filter applied to a dataset (***To be complete***)
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Printing out a filter should give us something at least mildly sensible.
setMethod("show",signature("filter"),function(object) 
	cat(paste("A filter named '",object@filterId,"'\n",sep="")))

## Most filters define only a single population
## setMethod("length","filter",function(x) 1)
setMethod("identifier", signature="filter",
          definition=function(object) object@filterId)


## ==========================================================================
## The %subset% operator constructs a filter that takes the subset of the LHS
## filter in the RHS filter result. 
## For many cases this is equivalent to an intersection filter.
## --------------------------------------------------------------------------
setMethod("%subset%",signature("filter","filter"),
          function(e1,e2) {
            new("subsetFilter",parameters=unique(c(e1@parameters,e2@parameters)),
                filters=list(e1,e2),filterId=paste(e1@filterId,"in",e2@filterId))
})

## ==========================================================================
## The | operator returns an object representing the union of two filters.
## A simple optimization would be to linearize the union of a filter and
## another union filter.
## --------------------------------------------------------------------------
setMethod("|",signature("filter","filter"),
          function(e1,e2) {
            new("unionFilter",parameters=unique(c(e1@parameters,e2@parameters)),
                filters=list(e1,e2),filterId=paste(e1@filterId,"or",e2@filterId))
          })


## ==========================================================================
## The & operator returns an object representing the intersection of two
## filters. This is somewhat different from the %subset% operator because
## some filters depend on the  data and would return different results
## when applied to the full dataset.
## --------------------------------------------------------------------------
setMethod("&",signature("filter","filter"),
          function(e1,e2) {
            new("intersectFilter",parameters=unique(c(e1@parameters,e2@parameters)),
                filters=list(e1,e2),filterId=paste(e1@filterId,"and",e2@filterId))
          })



## ==========================================================================
## The ! operator returns an object that simply takes the complement of the
## filter it returns.
## --------------------------------------------------------------------------
setMethod("!",signature("filter"),function(e1) 
	new("complementFilter",parameters=e1@parameters,filters=list(e1),
            filterId=paste("not",e1@filterId)))


## ==========================================================================
## --------------------------------------------------------------------------
setMethod("%&%",signature("ANY","ANY"),function(e1,e2) e1 %subset% e2)
