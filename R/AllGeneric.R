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
## Generic for applyFilter
## ---------------------------------------------------------------------------
setGeneric("applyTransformation", function(transformation,flowObject) standardGeneric("applyTransformation"))
## ---------------------------------------------------------------------------


## ===========================================================================
## Generic for applyFilter
## ---------------------------------------------------------------------------
setGeneric("applyFilter", function(filter,flowObject,parent) standardGeneric("applyFilter"))
## ---------------------------------------------------------------------------

## ===========================================================================
## Generic for applyFilter
## ---------------------------------------------------------------------------
setGeneric("basicplot", function(x,y,...) standardGeneric("basicplot"))
## ---------------------------------------------------------------------------

## ===========================================================================
## Generic for spillover
## ---------------------------------------------------------------------------
setGeneric("spillover",function(x,...) standardGeneric("spillover"))
## ---------------------------------------------------------------------------

## ===========================================================================
## Generic for applyCompensation
## ---------------------------------------------------------------------------
setGeneric("applyCompensation", function(compensation,flowObject) standardGeneric("applyCompensation"))
## ---------------------------------------------------------------------------
