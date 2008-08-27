## ==========================================================================
## Add or Replacement method to complete or replace the results
## coming out of a filtering operation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## simple accessor for the whole slot
setMethod("filterDetails",
          signature("filterResult", "missing"),
          function(result, filterId)
      {
          res <- result@filterDetails
          if(length(res)==1)
              res <- res[[1]]
          return(res)
      })

## access only a single filter of the filterDetails list
setMethod("filterDetails",
          signature("filterResult", "ANY"),
          function(result,filterId)
      {
          result@filterDetails[[filterId]]
      })

## replace a single filter in the list
setReplaceMethod("filterDetails",
                 signature("filterResult", "character", value="ANY"),
                 function(result,filterId,...,value)
             {
                 result@filterDetails[[filterId]] = value
                 result
             })

## replace the filter, this actually calls summarizeFilter
setReplaceMethod("filterDetails",
                 signature("filterResult", "character", value="filter"),
                 function(result, filterId, ..., value)
             {
                 filterDetails(result, filterId) <-
                     summarizeFilter(result,value)
                 result
             })

## For setOperationFilters we need to strip things from the attributes
setReplaceMethod("filterDetails",
                 signature("filterResult", "character",
                           value="setOperationFilter"),
                 function(result, filterId, ..., value)
             {
                 details <- attr(result@subSet,'filterDetails')
                 for(i in names(details)) {
                     filterDetails(result,i) <- details[[i]]
                 }
                 ##Record ourselves for posterity
                 filterDetails(result,filterId) <-
                     summarizeFilter(result,value)
                 result
             })



## ==========================================================================
## accessor method for parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature("filterResult"),
          function(object)
          filterDetails(object)$parameters)



## ==========================================================================
## Summarize a filtering operation
## This information will go in the filterDetails slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## By default, we just add one list item, which is the filter
setMethod("summarizeFilter",
          signature("filterResult", "filter"),
          function(result, filter)
      {
          list(filter=filter)
      })

## Add parameters to the filterDetails list if we have them
setMethod("summarizeFilter",
          signature("filterResult", "parameterFilter"),
          function(result,filter)
      {
          ret <- callNextMethod()
          ret$parameters <- parameters(filter)
          ret
      })


## summarize a filer
setMethod("summary",
          signature("filterResult"),
          function(object, ...)
      {
          summary(filterDetails(object, object@filterId)$filter,
                  object, ...)
      })





## ==========================================================================
## Allows us to compare filterResults and flowFrames. This lets us check
## (and warn or stop) that a particular flowFrame generated a filterResult
## allowing us to use it for further processing.
## --------------------------------------------------------------------------
setMethod("==",
          signature("flowFrame", "filterResult"),
          function(e1,e2)
      {
          i1 <- identifier(e1)
          i2 <- e2@frameId
          (length(i1) == 0 || length(i2) == 0 || i1 == i2)
      })

## Does S4 do this for us automagically? I don't know!
setMethod("==",
          signature("filterResult", "flowFrame"),
          function(e1, e2) e2==e1)



## ==========================================================================
## identifier method for filterResult
## --------------------------------------------------------------------------
setMethod("identifier", signature="filterResult",
          definition=function(object) object@filterId)




## ==========================================================================
## subset a filterResult by filterId, this only makes sense for
## multipleFilterResults, here we return everything...
## --------------------------------------------------------------------------
setMethod("[[",
          signature("filterResult"),
          function(x, i, j, drop=FALSE)
      {
          if((is.character(i) && i != x@filterId) || (as.numeric(i) > 1))
              stop("filter index out of bounds")
          x
      })


setMethod("show",
          signature("filterResult"),
          function(object) 
          cat(paste("A filterResult produced by the filter named '",
                    object@filterId, "'\n", sep="")))

