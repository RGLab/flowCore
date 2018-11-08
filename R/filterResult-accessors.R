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

#' Obtain details about a filter operation
#' 
#' A filtering operation captures details about its metadata and stores it in a
#' \code{filterDetails} slot in a \code{\linkS4class{filterResult}} object 
#' that is accessed using the \code{filterDetails} method. Each set of metadata 
#' is indexed by the \code{filterId} of the filter allowing for all the metadata 
#' in a complex filtering operation to be recovered after the final filtering.
#' 
#' 
#' @name filterDetails-methods
#' @aliases filterDetails-methods filterDetails
#' filterDetails,filterResult,missing-method
#' filterDetails,filterResult,ANY-method filterDetails<-
#' filterDetails<-,filterResult,character,setOperationFilter-method
#' filterDetails<-,filterResult,character,filter-method
#' filterDetails<-,filterResult,character,ANY-method
#' @usage NULL
#' @docType methods
#' 
#' @section Methods: 
#' \describe{
#' 
#' \item{filterDetails(result = "filterResult", filterId = "missing")}{When no particular
#' \code{filterId} is specified all the details are returned}
#' 
#' \item{filterDetails(result = "filterResult", filterId = "ANY")}{You can also obtain a
#' particular subset of details} 
#' }
#' 
#' @author B. Ellis, P.D. Haaland and N. LeMeur
#' @keywords methods
#' @export
setMethod("filterDetails",
          signature=signature(result="filterResult",
                              filterId="missing"),
          definition=function(result, filterId)
      {
          return(result@filterDetails)
      })

## Access only a single filter of the filterDetails list either by name
## or by index. This is only useful for manyFilterResults.
#' @export
setMethod("filterDetails",
          signature=signature(result="filterResult",
                              filterId="ANY"),
          definition=function(result, filterId)
      {
          result@filterDetails[[filterId]]
      })

## Replace a single filterDetails item in the list.
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
setMethod("==",
          signature=signature(e1="filterResult",
                              e2="flowFrame"),
          definition=function(e1, e2) e2==e1)




## ==========================================================================
## Subset a filterResult by filterId, this only makes real sense for
## multipleFilterResults, in the fallback option here we return everything...
## --------------------------------------------------------------------------
#' @export
setMethod("[[",
          signature=signature(x="filterResult"),
          definition=function(x, i, j, drop=FALSE)
      {
          if((is.character(i) && i != x@filterId) || (as.numeric(i) > 1))
              stop("filter index out of bounds")
          x
      })




