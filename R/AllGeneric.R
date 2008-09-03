## =========================================================================##
## =========================================================================##
##                       Generic methods definition                         ##
## =========================================================================##
## =========================================================================##


## ===========================================================================
## Generics for our operators
## %in% asks a question whereas %subset% forms a set.
## (default for %in% in base) 
## ---------------------------------------------------------------------------
setGeneric("%in%",function(x,table) standardGeneric("%in%"))
setGeneric("%subset%",function(e1,e2) standardGeneric("%subset%"))
setGeneric("%&%",function(e1,e2) standardGeneric("%&%"))
setGeneric("%on%",function(e1,e2) standardGeneric("%on%"))


## ===========================================================================
## Generic for sorting (default in base)
## ---------------------------------------------------------------------------
setGeneric("sort",function(x,decreasing=FALSE,...) standardGeneric("sort"))


## ===========================================================================
## Generic for range should be present via the 'Summary' group generic
## ---------------------------------------------------------------------------


## ===========================================================================
## Generic for ncol and nrow (defaults in base)
## ---------------------------------------------------------------------------
setGeneric("ncol", function(x) standardGeneric("ncol"))
setGeneric("nrow", function(x) standardGeneric("nrow"))


## ===========================================================================
## Generic for head and tail (already S3 in utils)
## ---------------------------------------------------------------------------
setGeneric("head", function(x, ...) standardGeneric("head"))
setGeneric("tail", function(x, ...) standardGeneric("tail"))


## ===========================================================================
## Generic for keyword
## ---------------------------------------------------------------------------
setGeneric("keyword",function(object,keyword) standardGeneric("keyword"))
setGeneric("keyword<-",function(object,value) standardGeneric("keyword<-"))
## setGeneric("description",function(object, hideInternal=FALSE)
##            standardGeneric("description"))
## setGeneric("description<-",function(object, value)
##            standardGeneric("description<-"))



## ===========================================================================
## Generic to access parameters in a flowframe
## ---------------------------------------------------------------------------
setGeneric("parameters",function(object) standardGeneric("parameters"))
setGeneric("parameters<-",
           function(object, value) standardGeneric("parameters<-"))


## ===========================================================================
## Generic to access colnames (defaults in base)
## ---------------------------------------------------------------------------
setGeneric("colnames",function(x, do.NULL = TRUE, prefix = "col") standardGeneric("colnames"))
setGeneric("colnames<-",
           function(x, value) standardGeneric("colnames<-"))


## ===========================================================================
## Generic for identifier
## ---------------------------------------------------------------------------
setGeneric("identifier",function(object) standardGeneric("identifier"))
setGeneric("identifier<-",function(object,value)
           standardGeneric("identifier<-"))


## ===========================================================================
## Generic for transformation (already S3 in base)
## ---------------------------------------------------------------------------
setGeneric("transform")


## ===========================================================================
## Generics for split and Subset (split is already S3 in base)
## ---------------------------------------------------------------------------
setGeneric("split")
setGeneric("Subset",function(x,subset,...) standardGeneric("Subset"))


## ===========================================================================
## Generic for filter (masking function in stats)
## ---------------------------------------------------------------------------
setGeneric("filter", function(x, filter,...) standardGeneric("filter"),
           useAsDefault=FALSE)


## ===========================================================================
## Generic for filterReference
## ---------------------------------------------------------------------------
setGeneric("filterReference",function(from,name)
           standardGeneric("filterReference"))


## ===========================================================================
## Generics for accessing and replacing certain aspects of the
## filtering process.
## ---------------------------------------------------------------------------
setGeneric("filterDetails<-", function(result, filterId,...,value)
           standardGeneric("filterDetails<-"))
setGeneric("filterDetails", function(result, filterId,...)
           standardGeneric("filterDetails"))


## ===========================================================================
## Used to summarize the operation of a filter on a frame. Used for the
## implementation of filters, not for the use of filters.
## ---------------------------------------------------------------------------
setGeneric("summarizeFilter", function(result,filter)
           standardGeneric("summarizeFilter"))


## ===========================================================================
## Used to summarize the operation of a filter on a frame (already S3 in base)
## ---------------------------------------------------------------------------
setGeneric("summary", function(object,...) standardGeneric("summary"))


## ===========================================================================
## Generics for spillover and compensate
## ---------------------------------------------------------------------------
setGeneric("spillover", function(x,...) standardGeneric("spillover"))
setGeneric("compensate", function(x,spillover, inv=TRUE, ...)
           standardGeneric("compensate"))


## ===========================================================================
## Generics for apply-like methods
## ---------------------------------------------------------------------------
setGeneric("fsApply",function(x,FUN,...,simplify=TRUE,use.exprs=FALSE)
           standardGeneric("fsApply"))
setGeneric("each_col",function(x,FUN,...) standardGeneric("each_col"))
setGeneric("each_row",function(x,FUN,...) standardGeneric("each_row"))



## ===========================================================================
## Generics for toTable for output in table-like structure
## ---------------------------------------------------------------------------
setGeneric("toTable",function(x,...) standardGeneric("toTable"))
