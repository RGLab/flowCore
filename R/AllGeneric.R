##  Generic methods definition - last changes PDH November 2, 2006
##      changed subSet to parent in applyFilter



## ===========================================================================
## These are already defined as generic functions:
## ---------------------------------------------------------------------------
if(!isGeneric("colnames<-"))
  setGeneric("colnames<-",    function(x, value) standardGeneric("colnames<-"))

if(!isGeneric("colnames"))
  setGeneric("colnames",      function(x, do.NULL=TRUE, prefix="col") standardGeneric("colnames"))

if(!isGeneric("plot"))
  setGeneric("plot",          useAsDefault=plot)

if(!isGeneric("hist"))
  setGeneric("hist",          useAsDefault=hist)

if(!isGeneric("lines"))
  setGeneric("lines",          useAsDefault=lines)

#if(!isGeneric("%in%"))

## %in% asks a question whereas %subset% forms a set.
setGeneric("%in%",function(x,table) standardGeneric("%in%"))
setGeneric("%subset%",function(e1,e2) standardGeneric("%subset%"))
setGeneric("%&%",function(e1,e2) standardGeneric("%&%"))


setGeneric("keyword",function(object,keyword) standardGeneric("keyword"))


## ---------------------------------------------------------------------------

## ===========================================================================
## Generic for nrow
## ---------------------------------------------------------------------------
if(!isGeneric("nrow"))
  setGeneric("nrow", function(x) standardGeneric("nrow"))
## ---------------------------------------------------------------------------

## ===========================================================================
## Generic for nrow
## ---------------------------------------------------------------------------
if(!isGeneric("ncol"))
  setGeneric("ncol", function(x) standardGeneric("ncol"))
## ---------------------------------------------------------------------------

setGeneric("identifier",function(object) standardGeneric("identifier"))


## ===========================================================================
## Generic for transformation
## ---------------------------------------------------------------------------
  setGeneric("transform")
## ---------------------------------------------------------------------------

## ===========================================================================
## Generic for split and subset
## ---------------------------------------------------------------------------
setGeneric("split")
setGeneric("Subset",function(x,subset,...) standardGeneric("Subset"))
## ---------------------------------------------------------------------------


## ===========================================================================
## Generic for filter
## ---------------------------------------------------------------------------
setGeneric("filter", function(x, filter,...) standardGeneric("filter"),useAsDefault=FALSE)

## For accessing certain aspects of the filtering process.
setGeneric("filterDetails<-",function(result,filterId,...,value)  standardGeneric("filterDetails<-"))
setGeneric("filterDetails",function(result,filterId,...) standardGeneric("filterDetails"))

## Used to summarize the operation of a filter on a frame. Used for the implementation
## of filters, not for the use of filters.
setGeneric("summarizeFilter",function(result,filter) standardGeneric("summarizeFilter"))


## ---------------------------------------------------------------------------


## ===========================================================================
## Generic for spillover
## ---------------------------------------------------------------------------
setGeneric("spillover",function(x,...) standardGeneric("spillover"))
setGeneric("compensate",function(x,spillover) standardGeneric("compensate"))
## ---------------------------------------------------------------------------

setGeneric("fsApply",function(x,FUN,...) standardGeneric("fsApply"))
