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
setGeneric("%in%",function(x,table) standardGeneric("%in%"))

setGeneric("%subset%",function(e1,e2) standardGeneric("%subset%"))

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


## ===========================================================================
## Generic for applyTransformation
## ---------------------------------------------------------------------------
  setGeneric("transform")
## ---------------------------------------------------------------------------


## ===========================================================================
## Generic for filter
## ---------------------------------------------------------------------------
setGeneric("filter", function(flowObject, filter) standardGeneric("filter"),useAsDefault=FALSE)
setGeneric("filterDetails",function(flowObject,filter,result) standardGeneric("filterDetails"))
## ---------------------------------------------------------------------------


## ===========================================================================
## Generic for spillover
## ---------------------------------------------------------------------------
setGeneric("spillover",function(x,...) standardGeneric("spillover"))
setGeneric("compensate",function(x,spillover) standardGeneric("compensate"))
## ---------------------------------------------------------------------------
