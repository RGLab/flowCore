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



## ===========================================================================
## Generic to access parameters in a flowframe
## ---------------------------------------------------------------------------
setGeneric("parameters", function(object, ...) standardGeneric("parameters"))
setGeneric("parameters<-",
           function(object, value) standardGeneric("parameters<-"))



## ===========================================================================
## Generic to access colnames (defaults in base)
## ---------------------------------------------------------------------------
setGeneric("colnames",function(x, do.NULL = TRUE, prefix = "col")
           standardGeneric("colnames"))
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
setGeneric("get",function(x, pos=-1, envir=as.environment(pos),
                          mode="any", inherits=TRUE)
           standardGeneric("get"))
setGeneric("mget", function (x, envir, mode="any",
                             ifnotfound=list(function(x)
                                             stop(paste("value for '", 
                                                        x, "' not found",
                                                        sep = ""),
                                                  call.=FALSE)),
                             inherits=FALSE)
           standardGeneric("mget"))
setGeneric("assign", function(x, value, pos=-1, envir=as.environment(pos),
                              inherits=FALSE, immediate=TRUE)
           standardGeneric("assign"))



## ===========================================================================
## Generics to list the content of a workFlow object
## ---------------------------------------------------------------------------
setGeneric("ls", function(name, pos=-1, envir=as.environment(pos),
                          all.names=FALSE, pattern)
            standardGeneric("ls"))



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
## Generic function eval for all the transformations
## ---------------------------------------------------------------------------
setGeneric("eval",
            def=function(expr,
                        envir = parent.frame(),
                        enclos = if (is.list(envir) ||is.pairlist(envir))    
                                      parent.frame() 
                                 else 
                                      baseenv()) 
            standardGeneric("eval"),
            useAsDefault=function(expr,envir=new.env(),enclos)
 	    {
           	stop("'eval' not defined on class '",class(object), "'")
            }
 	  )




## ===========================================================================
## Generic to apply a normalization object on a flowSet
## ---------------------------------------------------------------------------
setGeneric("normalize", function(data, x)
           standardGeneric("normalize"))
