## =========================================================================##
## =========================================================================##
##                       Generic methods definitions                        ##
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
## Generic for printing details of an object
## ---------------------------------------------------------------------------
setGeneric("print",function(x,...) standardGeneric("print"))



## ===========================================================================
## Generic for range should be present via the 'Summary' group generic
## ---------------------------------------------------------------------------



## ===========================================================================
## Generic for head and tail (already S3 in utils)
## ---------------------------------------------------------------------------
setGeneric("head", function(x, ...) standardGeneric("head"))
setGeneric("tail", function(x, ...) standardGeneric("tail"))



## ===========================================================================
## Generic for keyword
## ---------------------------------------------------------------------------
setGeneric("keyword",function(object,keyword, ...) standardGeneric("keyword"))
setGeneric("keyword<-",function(object,value) standardGeneric("keyword<-"))



## ===========================================================================
## Generic to access parameters in a flowframe
## ---------------------------------------------------------------------------
setGeneric("parameters", function(object, ...) standardGeneric("parameters"))
setGeneric("parameters<-",
           function(object, value) standardGeneric("parameters<-"))



## ===========================================================================
## Generic to access colnames (defaults in base)
## ---------------------------------------------------------------------------
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
## Generic for filter (using stats::filter)
## ---------------------------------------------------------------------------
setGeneric("filter")



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
## Generics for spillover and compensate
## ---------------------------------------------------------------------------
setGeneric("spillover", function(x,...) standardGeneric("spillover"))
setGeneric("compensate", 
            function(x, spillover, ...)
            standardGeneric("compensate"))


## ===========================================================================
## Generics for apply-like methods
## ---------------------------------------------------------------------------
setGeneric("fsApply",function(x,FUN,...,simplify=TRUE,use.exprs=FALSE)
           standardGeneric("fsApply"))
setGeneric("each_col",function(x,FUN,...) standardGeneric("each_col"))
setGeneric("each_row",function(x,FUN,...) standardGeneric("each_row"))



## ===========================================================================
## Generic for toTable for output in table-like structure
## ---------------------------------------------------------------------------
setGeneric("toTable",function(x,...) standardGeneric("toTable"))



## ===========================================================================
## Generic to test for a NULL reference
## setMethod("is.null",
##           signature=signature("fcReference"),
##           definition=function(f) is(f, "fcNullReference"))
## Doesn't work in this R version. We are not allowed to define methods
## on all primitives
## ---------------------------------------------------------------------------
setGeneric("isNull",function(f) standardGeneric("isNull"))



## ===========================================================================
## Generics for parentView and gate slot accessors (of actionItems)
## ---------------------------------------------------------------------------
setGeneric("parent",function(object) standardGeneric("parent"))
setGeneric("gate",function(object) standardGeneric("gate"))



## ===========================================================================
## Generics for action slot accessors (of views)
## ---------------------------------------------------------------------------
setGeneric("action",function(object) standardGeneric("action"))



## ===========================================================================
## Generics to add an actionItem to a workflow
## ---------------------------------------------------------------------------
setGeneric("add",function(wf, action, ...) standardGeneric("add"))



## ===========================================================================
## Generics to assign values into and get values from a workflow environment
## ---------------------------------------------------------------------------
setGeneric("assign", function(x, value, pos=-1, envir=as.environment(pos),
                              inherits=FALSE, immediate=TRUE)
           standardGeneric("assign"))



## ===========================================================================
## Generics to list the content of a workFlow object(using base::ls)
## ---------------------------------------------------------------------------
setGeneric("ls")



## ===========================================================================
## Generics to remove object from a workFlow. We pass along a substituted
## version of symbol in order to be able to remove references as well.
## ---------------------------------------------------------------------------
setGeneric("Rm", function(symbol, envir, subSymbol, ...)
       {
           subSymbol <- as.character(substitute(symbol))
           standardGeneric("Rm")
       })



## ===========================================================================
## Generic to get to the data linked to a view or actionItem
## ---------------------------------------------------------------------------
setGeneric("Data", function(object)
           standardGeneric("Data"))



## ===========================================================================
## Generic to access the alias table of a workFlow
## ---------------------------------------------------------------------------
setGeneric("alias", function(object, ...)
           standardGeneric("alias"))



## ===========================================================================
## Generic to access the journal of a workFlow
## ---------------------------------------------------------------------------
setGeneric("journal", function(object, ...)
           standardGeneric("journal"))



## ===========================================================================
## Generic to access the tree within a workFlow
## ---------------------------------------------------------------------------
setGeneric("tree", function(object, ...)
           standardGeneric("tree"))



## ===========================================================================
## Generic to return all views of the workFlow object
## ---------------------------------------------------------------------------
setGeneric("views", function(x, ...)
           standardGeneric("views"))



## ===========================================================================
## Generic to return all actionItems of the workFlow object
## ---------------------------------------------------------------------------
setGeneric("actions", function(x, ...)
           standardGeneric("actions"))



## ===========================================================================
## Generic to apply a normalization object on a flowSet
## ---------------------------------------------------------------------------
setGeneric("normalize", function(data, x,...)
           standardGeneric("normalize"))

## ===========================================================================
## Generic to extract index sorted data
## ---------------------------------------------------------------------------
setGeneric("getIndexSort",def=function(x){
  standardGeneric("getIndexSort")
})

