## ==========================================================================
## summarizeFilter methods collect all necessary information from a
## filtering operation and store them in the filterDetails list. In the
## most basic case, this is simply the filter object. parameterFilters
## should all add a parameters item. The return value of these methods
## has to be a list.
## ==========================================================================


#' Methods for function summarizeFilter
#' 
#' 
#' Internal methods to populate the \code{filterDetails} slot of a
#' \code{\link{filterResult}} object.
#' 
#' 
#' @name summarizeFilter-methods
#' @aliases summarizeFilter-methods summarizeFilter
#' summarizeFilter,filterResult,filter-method
#' summarizeFilter,filterResult,filterReference-method
#' summarizeFilter,filterResult,parameterFilter-method
#' summarizeFilter,filterResult,subsetFilter-method
#' summarizeFilter,logicalFilterResult,norm2Filter-method
#' summarizeFilter,logicalFilterResult,parameterFilter-method
#' summarizeFilter,multipleFilterResult,parameterFilter-method
#' @usage summarizeFilter(result, filter)
#' @docType methods
#' @param result A \code{\linkS4class{filterResult}} (or one of its derived classes)
#' representing the result of a filtering operation in whose 
#' \code{filterDetails} slot the information will be stored.
#' @param filter The corresponding \code{\linkS4class{filter}} (or one of its
#' derived classes). 
#' @section Methods: 
#' \describe{
#' 
#' \item{summarizeFilter(result = "filterResult", filter = "filter")}{ \code{summarizeFilter}
#' methods are called during the process of filtering. Their output is a list,
#' and it can be arbitrary data that should be stored along with the results of
#' a filtering operation. }
#' 
#' \item{summarizeFilter(result = "filterResult", filter = "filterReference")}{ see above }
#' 
#' \item{summarizeFilter(result = "filterResult", filter = "parameterFilter")}{ see above }
#' 
#' \item{summarizeFilter(result = "filterResult", filter = "subsetFilter")}{ see above }
#' 
#' \item{summarizeFilter(result = "logicalFilterResult", filter = "norm2Filter")}{ see above }
#' 
#' \item{summarizeFilter(result = "logicalFilterResult", filter = "parameterFilter")}{ see above
#' }
#' 
#' \item{summarizeFilter(result = "multipleFilterResult", filter = "parameterFilter")}{ see
#' above } }
#' @keywords methods


## ==========================================================================
## By default, we just add a single list item, which is the filter
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature(result="filterResult",
                    filter="filter"),
          definition=function(result, filter)
      {
          list(filter=filter)
      })



## ==========================================================================
## Add parameters to the filterDetails list if we have them.
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature=signature(result="filterResult",
                              filter="parameterFilter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$parameters <- parameters(filter)
          ret
      })



## ==========================================================================
## For filterReferences, we first resolve to the concrete filter and
## summarize that
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature=signature(result="filterResult",
                              filter="filterReference"),
          definition=function(result, filter)
          summarizeFilter(result, as(filter, "concreteFilter")))



## ==========================================================================
## The default for multipleFilterResults is to take the levels in the subSet
## slot as the population item.
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature=signature(result="multipleFilterResult",
                              filter="parameterFilter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$populations <- levels(result@subSet)
          return(ret)
      })



## ==========================================================================
## The default for logicalFilterResults is to create useful population names
## from the filter name
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature=signature(result="logicalFilterResult",
                              filter="parameterFilter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$populations <- paste(identifier(filter), c("+", "-"), sep="")
          return(ret)
      })


## ==========================================================================
## For a norm2Filter we want to strip things from the attributes in the
## subSet slot, i.e., the details about the fitted bivariate normal
## distribution
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature=signature(result="logicalFilterResult",
                              filter="norm2Filter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$cov <- attr(result@subSet,'cov')
          ret$center <- attr(result@subSet,'center')
          ret$radius <- attr(result@subSet,'radius')
          ret$parameters <- parameters(filter)
          return(ret)
      })



## ==========================================================================
## For a subsetFilter we need to grab things from the attributes
## ---------------------------------------------------------------------------
#' @export
setMethod("summarizeFilter",
          signature=signature(result="filterResult",
                              filter="subsetFilter"),
          definition=function(result,filter)
      {
          ret <- callNextMethod()
          ret$subsetCount <- attr(result@subSet,'subsetCount')
          additionalDetails <- attr(result@subSet,'filterDetails')[[2]]
          sel <- !names(additionalDetails) %in% c("filter", "parameters")
          ret <- c(ret, additionalDetails[sel])
          return(ret)
      })
