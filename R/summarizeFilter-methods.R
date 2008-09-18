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
          signature("filterResult", "filter"),
          definition=function(result, filter)
      {
          list(filter=filter)
      })



## ==========================================================================
## Add parameters to the filterDetails list if we have them.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature("filterResult", "parameterFilter"),
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
          signature=signature("filterResult", "filterReference"),
          definition=function(result, filter)
          summarizeFilter(result, as(filter, "concreteFilter")))



## ==========================================================================
## The default for multipleFilterResults is to take the levels in the subSet
## slot as the population item.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature("multipleFilterResult", "parameterFilter"),
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
          signature=signature("logicalFilterResult", "parameterFilter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$populations <- paste(identifier(filter), c("+", "-"), sep="")
          return(ret)
      })



## ==========================================================================
## For curv1Filters we want to strip the attributes from the subSet slot,
## i.e., the boundaries of the high density regions and the fsObj obects.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature("multipleFilterResult","curv1Filter"),
          definition=function(result,filter)
      {
          ret <- callNextMethod()
          ret$boundaries <- attr(result@subSet, "boundaries")
          ret$fsObj <- attr(result@subSet, "fSObj")
          return(ret)
      })



## ==========================================================================
## For curv2Filters we want to strip the attributes from the subSet slot,
## i.e., the boundaries of the high density regions and the fsObj obects.
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature("multipleFilterResult","curv2Filter"),
          definition=function(result,filter)
      {
          ret <- callNextMethod()
          ret$polygons <- attr(result@subSet, "polygon")
	  ret$fsObj <- attr(result@subSet, "fSObj")
          return(ret)
      })



## ==========================================================================
## For a norm2Filter we want to strip things from the attributes in the
## subSet slot, i.e., the details about the fitted bivariate normal
## distribution
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature("logicalFilterResult", "norm2Filter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          ret$cov <- attr(result@subSet,'cov')
          ret$center <- attr(result@subSet,'center')
          ret$radius <- attr(result@subSet,'radius')
          return(ret)
      })



## ==========================================================================
## For a subsetFilter we need to grab things from the attributes
## ---------------------------------------------------------------------------
setMethod("summarizeFilter",
          signature=signature("filterResult","subsetFilter"),
          definition=function(result,filter)
      {
          ret <- callNextMethod()
          ret$subsetCount <- attr(result@subSet,'subsetCount')
          return(ret)
      })
