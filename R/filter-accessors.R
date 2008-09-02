## ==========================================================================
## This implementation of filter does not require a per-filter implementation
## anything filter specific has been moved to either filterDetails or the %in%
## implementation. A filter operation is now about assembling those two calls 
## AND associating a filterResult with a particular flowFrame if possible.
## ==========================================================================



## ==========================================================================
## accessor method for slot parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature("filter"),
          definition=function(object) character(0))

setMethod("parameters",
          signature=signature("parameterFilter"),
          definition=function(object) object@parameters)

setMethod("parameters",
          signature=signature("setOperationFilter"),
          definition=function(object)
          unique(unlist(lapply(object@filters, parameters))))



## ==========================================================================
## summary method for filters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("summary",
          signature=signature("filter"),
          definition=function(object, result, ...)
      {
          if(missing(result))
              stop("Only resolved filters may be summarized.\nTry something ",
                   "like 'filter(myFlowFrame, myFilter)'\nbefore calling ",
                   "'summary'.", call.=FALSE)
          else
              l <- as(result, "logical")
          true <- sum(l)
          count <- length(l)
          new("filterSummary", name=identifier(object), true=true,
              count=count, p=true/count)
      })



## ==========================================================================
## Printing out a filter should give us something at least mildly sensible.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
          signature=signature("filter"),
          definition=function(object) 
          cat(paste("A filter named '", object@filterId, "'\n", sep="")))


## Most filters define only a single population
## setMethod("length","filter",function(x) 1)



## ==========================================================================
## We need to be able to get an ID for a filter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("identifier",
          signature="filter",
          definition=function(object) object@filterId)

setReplaceMethod("identifier",
                 signature=signature("filter","character"),
                 definition=function(object,value)
             {
                 object@filterId = value
                 object
             })



## ==========================================================================
## The %subset% operator constructs a filter that takes the subset of the LHS
## filter in the RHS filter result. 
## For many cases this is equivalent to an intersection filter.
## --------------------------------------------------------------------------
setMethod("%subset%",
          signature=signature("filter","filter"),
          definition=function(e1, e2)
      {
          new("subsetFilter",
              filters=list(e1, e2), filterId=paste(identifier(e1),"in",
                                    identifier(e2)))
      })

setMethod("%subset%",
          signature=signature("list","filter"),
          definition=function(e1, e2) lapply(e1, "%subset%", e2=e2))

setMethod("%subset%",
          signature=signature("filterSet","filter"),
          definition=function(e1,e2)
      {
          ## Make a copy of the filterSet, preserving R semantics
          x <- as(as(e1, "list"), "filterSet")
          n <- names(e1)
          x[[""]] <- e2
          target <- as.symbol(identifier(e2))
          for(i in n){
              x[[""]] <- substitute(~ a %subset% b, list(a=as.symbol(i),
                                                         b=target))
          }
          x
      })



## ==========================================================================
## The | operator returns an object representing the union of two filters.
## A simple optimization would be to linearize the union of a filter and
## another union filter.
## --------------------------------------------------------------------------
setMethod("|",
          signature=signature("filter", "filter"),
          definition=function(e1, e2)
      {
          new("unionFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "or", identifier(e2)))
      })

setMethod("|",
          signature=signature("list","filter"),
          definition=function(e1, e2) lapply(e1, "|", e2=e2))

setMethod("|",
          signature=signature("filter", "list"),
          definition=function(e1, e2) lapply(e2, "|", e1=e1))



## ==========================================================================
## The & operator returns an object representing the intersection of two
## filters. This is somewhat different from the %subset% operator because
## some filters depend on the data and would return different results
## when applied to the full dataset.
## --------------------------------------------------------------------------
setMethod("&",
          signature=signature("filter", "filter"),
          definition=function(e1, e2)
      {
          new("intersectFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "and", identifier(e2)))
      })

setMethod("&",
          signature=signature("list", "filter"),
          definition=function(e1, e2) lapply(e1, "&", e2=e2))

setMethod("&",
          signature=signature("filter", "list"),
          definition=function(e1, e2) lapply(e2, "&", e1=e1))



## ==========================================================================
## The ! operator returns an object that simply takes the complement of the
## filter it returns.
## --------------------------------------------------------------------------
setMethod("!",
          signature=signature("filter"),
          definition=function(x)
      {
	new("complementFilter",filters=list(x),
            filterId=paste("not",identifier(x)))
})



## ==========================================================================
## --------------------------------------------------------------------------
setMethod("%&%",
          signature=signature("filter", "filter"),
          definition=function(e1, e2) e1 %subset% e2)

