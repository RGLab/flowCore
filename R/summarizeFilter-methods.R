## ==========================================================================
## summarizeFilter methods collect all necessary information from a
## filtering operation and store them in the filterDetails list. In the
## most basic case, this is simply the filter object. parameterFilters
## should all add a parameters item. The return value of these methods
## has to be a list.
## ==========================================================================




## ==========================================================================
## By default, we just add a single list item, which is the filter
## ---------------------------------------------------------------------------
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
setMethod("summarizeFilter",
          signature=signature(result="filterResult",
                              filter="filterReference"),
          definition=function(result, filter)
          summarizeFilter(result, as(filter, "concreteFilter")))



## ==========================================================================
## The default for multipleFilterResults is to take the levels in the subSet
## slot as the population item.
## ---------------------------------------------------------------------------
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
setMethod("summarizeFilter",
          signature=signature(result="logicalFilterResult",
                              filter="parameterFilter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$populations <- paste(identifier(filter), c("+", "-"), sep="")
          return(ret)
      })


#TODO: to remove
## ==========================================================================
## For curv1Filters we want to strip the attributes from the subSet slot,
## i.e., the boundaries of the high density regions and the fsObj obects.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature(result="multipleFilterResult",
                              filter="curv1Filter"),
          definition=function(result,filter)
      {
          ret <- callNextMethod()
          ret$boundaries <- attr(result@subSet, "boundaries")
          ret$fsObj <- attr(result@subSet, "fSObj")
          return(ret)
      })


#TODO: to remove
  
## ==========================================================================
## For curv2Filters we want to strip the attributes from the subSet slot,
## i.e., the boundaries of the high density regions and the fsObj obects.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature(result="multipleFilterResult",
                              filter="curv2Filter"),
          definition=function(result,filter)
      {
          ret <- callNextMethod()
          ret$polygons <- attr(result@subSet, "polygon")
	  ret$fsObj <- attr(result@subSet, "fSObj")
          return(ret)
      })


#TODO: to remove
  
## ==========================================================================
## For a norm2Filter we want to strip things from the attributes in the
## subSet slot, i.e., the details about the fitted bivariate normal
## distribution
## ---------------------------------------------------------------------------
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
