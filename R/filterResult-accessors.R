## ==========================================================================
## filterResults are the output of a filtering operation. 
## ==========================================================================






## ==========================================================================
## Add or replacement method to complete or replace the results coming out
## of a filtering operation. The filter replacement method (along with it
## the summarizeFilter method) gets called during each filtering operation,
## making sure that the filter item in the filterDetails slot contains all
## necessary information
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Simple accessor for the whole slot. Most filterResult classes except for
## manyFilterResults only contain a single filterDetails item, and those are
## returned directly.
setMethod("filterDetails",
          signature=signature(result="filterResult",
                              filterId="missing"),
          definition=function(result, filterId)
      {
          return(result@filterDetails)
      })

## Access only a single filter of the filterDetails list either by name
## or by index. This is only useful for manyFilterResults.
setMethod("filterDetails",
          signature=signature(result="filterResult",
                              filterId="ANY"),
          definition=function(result, filterId)
      {
          result@filterDetails[[filterId]]
      })

## Replace a single filterDetails item in the list. 
setReplaceMethod("filterDetails",
                 signature=signature(result="filterResult",
                                     filterId="character",
                 value="ANY"),
                 definition=function(result, filterId, ... ,value)
             {
                 result@filterDetails[[filterId]] <- value
                 result
             })

## Replace the filter item in a single filterDetails list. This is actually
## only a wrapper around summarizeFilter which has to be defined for each
## filter class separately unless the default behaviour of simply adding
## the filter is sufficient
setReplaceMethod("filterDetails",
                 signature=signature(result="filterResult",
                                     filterId="character",
                 value="filter"),
                 definition=function(result, filterId, ..., value)
             {
                 filterDetails(result, filterId) <-
                     summarizeFilter(result,value)
                 result
             })

## For setOperationFilters we need to strip the information for the
## individual filters in the manyFilterResult from the attributes
setReplaceMethod("filterDetails",
                 signature=signature(result="filterResult",
                                     filterId="character",
                 value="setOperationFilter"),
                 definition=function(result, filterId, ..., value)
             {
                 details <- attr(result@subSet, 'filterDetails')
                 for(i in names(details)) {
                     filterDetails(result, i) <- details[[i]]
                 }
                 filterDetails(result, filterId) <-
                     summarizeFilter(result, value)
                 result
             })



## ==========================================================================
## Allows us to compare filterResults and flowFrames. This lets us check
## (and warn or stop) that a particular flowFrame generated a filterResult
## allowing us to use it for further processing.
## --------------------------------------------------------------------------
setMethod("==",
          signature=signature(e1="flowFrame",
                              e2="filterResult"),
          definition=function(e1, e2)
      {
          i1 <- identifier(e1)
          i2 <- e2@frameId
          (length(i1) == 0 || length(i2) == 0 || i1 == i2)
      })

## Does S4 do this for us automagically? I don't know!
setMethod("==",
          signature=signature(e1="filterResult",
                              e2="flowFrame"),
          definition=function(e1, e2) e2==e1)




## ==========================================================================
## Subset a filterResult by filterId, this only makes real sense for
## multipleFilterResults, in the fallback option here we return everything...
## --------------------------------------------------------------------------
setMethod("[[",
          signature=signature(x="filterResult"),
          definition=function(x, i, j, drop=FALSE)
      {
          if((is.character(i) && i != x@filterId) || (as.numeric(i) > 1))
              stop("filter index out of bounds")
          x
      })




