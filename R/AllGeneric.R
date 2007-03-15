##  Generic methods definition - last changes NLM January 17, 2006
##      add parameters

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

setGeneric("draw",function(x,data,...) standardGeneric("draw"))

## ===========================================================================
## Generics for our operators
## %in% asks a question whereas %subset% forms a set.
## ---------------------------------------------------------------------------
setGeneric("%in%",function(x,table) standardGeneric("%in%"))
setGeneric("%subset%",function(e1,e2) standardGeneric("%subset%"))
setGeneric("%&%",function(e1,e2) standardGeneric("%&%"))
setGeneric("%on%",function(e1,e2) standardGeneric("%on%"))


## ===========================================================================
## Generic for keyword
## ---------------------------------------------------------------------------
setGeneric("keyword",function(object,keyword) standardGeneric("keyword"))


## ===========================================================================
## Generic to access parameters in a flowframe
## ---------------------------------------------------------------------------
setGeneric("parameters",function(object) standardGeneric("parameters"))


## ===========================================================================
## Generic for nrow
## ---------------------------------------------------------------------------
if(!isGeneric("nrow"))
  setGeneric("nrow", function(x) standardGeneric("nrow"))


## ===========================================================================
## Generic for ncol
## ---------------------------------------------------------------------------
if(!isGeneric("ncol"))
  setGeneric("ncol", function(x) standardGeneric("ncol"))


## ===========================================================================
## Generic for identifier
## ---------------------------------------------------------------------------
setGeneric("identifier",function(object) standardGeneric("identifier"))


## ===========================================================================
## Generic for transformation
## ---------------------------------------------------------------------------
  setGeneric("transform")


## ===========================================================================
## Generics for split and Subset
## ---------------------------------------------------------------------------
setGeneric("split")
setGeneric("Subset",function(x,subset,...) standardGeneric("Subset"))



## ===========================================================================
## Generic for filter
## ---------------------------------------------------------------------------
setGeneric("filter", function(x, filter,...) standardGeneric("filter"),useAsDefault=FALSE)


## ===========================================================================
## Generics for accessing certain aspects of the filtering process.
## ---------------------------------------------------------------------------
setGeneric("filterDetails<-",function(result,filterId,...,value)  standardGeneric("filterDetails<-"))
setGeneric("filterDetails",function(result,filterId,...) standardGeneric("filterDetails"))


## ===========================================================================
## Used to summarize the operation of a filter on a frame. Used for the implementation
## of filters, not for the use of filters.
## ---------------------------------------------------------------------------
setGeneric("summarizeFilter",function(result,filter) standardGeneric("summarizeFilter"))


## ===========================================================================
## Generics for spillover and compensate
## ---------------------------------------------------------------------------
setGeneric("spillover",function(x,...) standardGeneric("spillover"))
setGeneric("compensate",function(x,spillover) standardGeneric("compensate"))


## ===========================================================================
## Generics for apply-like methods
## ---------------------------------------------------------------------------
setGeneric("fsApply",function(x,FUN,...,simplify=TRUE,use.exprs=FALSE) standardGeneric("fsApply"))
setGeneric("each_col",function(x,FUN,...) standardGeneric("each_col"))
setGeneric("each_row",function(x,FUN,...) standardGeneric("each_row"))
