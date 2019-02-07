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
#' @export
setGeneric("%in%",function(x,table) standardGeneric("%in%"))
#' @export
setGeneric("%subset%",function(e1,e2) standardGeneric("%subset%"))
#' @export
setGeneric("%&%",function(e1,e2) standardGeneric("%&%"))
#' @export
setGeneric("%on%",function(e1,e2) standardGeneric("%on%"))


## ===========================================================================
## Generic for printing details of an object
## ---------------------------------------------------------------------------
#' @export
setGeneric("print",function(x,...) standardGeneric("print"))



## ===========================================================================
## Generic for range should be present via the 'Summary' group generic
## ---------------------------------------------------------------------------



## ===========================================================================
## Generic for head and tail (already S3 in utils)
## ---------------------------------------------------------------------------
#' @export
setGeneric("head", function(x, ...) standardGeneric("head"))
#' @export
setGeneric("tail", function(x, ...) standardGeneric("tail"))



## ===========================================================================
## Generic for keyword
## ---------------------------------------------------------------------------
#' @export
setGeneric("keyword",function(object,keyword, ...) standardGeneric("keyword"))
#' @export
setGeneric("keyword<-",function(object,value) standardGeneric("keyword<-"))



## ===========================================================================
## Generic to access parameters in a flowframe
## ---------------------------------------------------------------------------
#' @export
setGeneric("parameters", function(object, ...) standardGeneric("parameters"))
#' @export
setGeneric("parameters<-",
           function(object, value) standardGeneric("parameters<-"))



## ===========================================================================
## Generic to access colnames (defaults in base)
## ---------------------------------------------------------------------------
#' @export
setGeneric("colnames<-",
           function(x, value) standardGeneric("colnames<-"))



## ===========================================================================
## Generic for identifier
## ---------------------------------------------------------------------------
#' @export
setGeneric("identifier",function(object) standardGeneric("identifier"))
#' @export
setGeneric("identifier<-",function(object,value)
           standardGeneric("identifier<-"))



## ===========================================================================
## Generic for transformation (already S3 in base)
## ---------------------------------------------------------------------------
#' @export
setGeneric("transform")



## ===========================================================================
## Generics for split and Subset (split is already S3 in base)
## ---------------------------------------------------------------------------
#' @export
setGeneric("split")
#' @export
setGeneric("Subset",function(x,subset,...) standardGeneric("Subset"))



## ===========================================================================
## Generic for filter (using stats::filter)
## ---------------------------------------------------------------------------
#' @export
setGeneric("filter")



## ===========================================================================
## Generic for filterReference
## ---------------------------------------------------------------------------
#' @export
setGeneric("filterReference",function(from,name)
           standardGeneric("filterReference"))



## ===========================================================================
## Generics for accessing and replacing certain aspects of the
## filtering process.
## ---------------------------------------------------------------------------
#' @export
setGeneric("filterDetails<-", function(result, filterId,...,value)
           standardGeneric("filterDetails<-"))
#' @export
setGeneric("filterDetails", function(result, filterId,...)
           standardGeneric("filterDetails"))



## ===========================================================================
## Used to summarize the operation of a filter on a frame. Used for the
## implementation of filters, not for the use of filters.
## ---------------------------------------------------------------------------
#' @export
setGeneric("summarizeFilter", function(result,filter)
           standardGeneric("summarizeFilter"))



## ===========================================================================
## Generics for spillover and compensate and decompensate
## ---------------------------------------------------------------------------
#' @export
setGeneric("spillover", function(x,...) standardGeneric("spillover"))
#' @export
setGeneric("compensate", 
            function(x, spillover, ...)
            standardGeneric("compensate"))
#' @export
setGeneric("decompensate",
			  function(x, spillover, ...)
			  	standardGeneric("decompensate"))

## ================================
## Generic for spillover_match
## --------------------------------
#' @export
setGeneric("spillover_match", function(x,...) standardGeneric("spillover_match"))

## ===========================================================================
## Generics for apply-like methods
## ---------------------------------------------------------------------------
#' @export
setGeneric("fsApply",function(x,FUN,...,simplify=TRUE,use.exprs=FALSE)
           standardGeneric("fsApply"))

#' Methods to apply functions over flowFrame margins
#' 
#' 
#' Returns a vector or array of values obtained by applying a function to the
#' margins of a flowFrame. This is equivalent of running \code{\link{apply}} on
#' the output of \code{exprs(flowFrame)}.
#' 
#' 
#' @name each_col
#' @aliases each_row each_row-methods each_row,flowFrame-method each_col
#' each_col-methods each_col,flowFrame-method
#' @docType methods
#' 
#' @usage 
#' each_col(x, FUN, ...)
#' each_row(x, FUN, ...)
#' 
#' @param x Object of class \code{\linkS4class{flowFrame}}.
#' @param FUN the function to be applied. In the case of functions like '+',
#' '\%*\%', etc., the function name must be backquoted or quoted.
#' @param ... optional arguments to 'FUN'.
#' 
#' 
#' @author B. Ellis, N. LeMeur, F. Hahne
#' @seealso \code{\link{apply}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata", "0877408774.B08", package="flowCore"),
#' transformation="linearize")
#' each_col(samp, summary)
#' 
#' 
#' @export
setGeneric("each_col",function(x,FUN,...) standardGeneric("each_col"))
#' @export
setGeneric("each_row",function(x,FUN,...) standardGeneric("each_row"))



## ===========================================================================
## Generic for toTable for output in table-like structure
## ---------------------------------------------------------------------------
#' @export
setGeneric("toTable",function(x,...) standardGeneric("toTable"))



## ===========================================================================
## Generic to test for a NULL reference
## setMethod("is.null",
##           signature=signature("fcReference"),
##           definition=function(f) is(f, "fcNullReference"))
## Doesn't work in this R version. We are not allowed to define methods
## on all primitives
## ---------------------------------------------------------------------------
#' @export
setGeneric("isNull",function(f) standardGeneric("isNull"))



## ===========================================================================
## Generics for parentView and gate slot accessors (of actionItems)
## ---------------------------------------------------------------------------
#' @export
setGeneric("parent",function(object) standardGeneric("parent"))
#' @export
setGeneric("gate",function(object) standardGeneric("gate"))



## ===========================================================================
## Generics for action slot accessors (of views)
## ---------------------------------------------------------------------------
#' @export
setGeneric("action",function(object) standardGeneric("action"))



## ===========================================================================
## Generics to add an actionItem to a workflow
## ---------------------------------------------------------------------------
#' @export
setGeneric("add",function(wf, action, ...) standardGeneric("add"))



## ===========================================================================
## Generics to assign values into and get values from a workflow environment
## ---------------------------------------------------------------------------
#' @export
setGeneric("assign", function(x, value, pos=-1, envir=as.environment(pos),
                              inherits=FALSE, immediate=TRUE)
           standardGeneric("assign"))



## ===========================================================================
## Generics to list the content of a workFlow object(using base::ls)
## ---------------------------------------------------------------------------
#' @export
setGeneric("ls")



## ===========================================================================
## Generics to remove object from a workFlow. We pass along a substituted
## version of symbol in order to be able to remove references as well.
## ---------------------------------------------------------------------------
#' @export
setGeneric("Rm", function(symbol, envir, subSymbol, ...)
       {
           subSymbol <- as.character(substitute(symbol))
           standardGeneric("Rm")
       })



## ===========================================================================
## Generic to get to the data linked to a view or actionItem
## ---------------------------------------------------------------------------
#' @export
setGeneric("Data", function(object)
           standardGeneric("Data"))



## ===========================================================================
## Generic to access the alias table of a workFlow
## ---------------------------------------------------------------------------
#' @export
setGeneric("alias", function(object, ...)
           standardGeneric("alias"))



## ===========================================================================
## Generic to access the journal of a workFlow
## ---------------------------------------------------------------------------
#' @export
setGeneric("journal", function(object, ...)
           standardGeneric("journal"))



## ===========================================================================
## Generic to access the tree within a workFlow
## ---------------------------------------------------------------------------
#' @export
setGeneric("tree", function(object, ...)
           standardGeneric("tree"))



## ===========================================================================
## Generic to return all views of the workFlow object
## ---------------------------------------------------------------------------
#' @export
setGeneric("views", function(x, ...)
           standardGeneric("views"))



## ===========================================================================
## Generic to return all actionItems of the workFlow object
## ---------------------------------------------------------------------------
#' @export
setGeneric("actions", function(x, ...)
           standardGeneric("actions"))



## ===========================================================================
## Generic to apply a normalization object on a flowSet
## ---------------------------------------------------------------------------
#' @export
setGeneric("normalize", function(data, x,...)
           standardGeneric("normalize"))

## ===========================================================================
## Generic to extract index sorted data
## ---------------------------------------------------------------------------
#' @export
setGeneric("getIndexSort",def=function(x){
  standardGeneric("getIndexSort")
})

## ===========================================================================
## Generics for modifying gates
## ---------------------------------------------------------------------------
#' @export
transform_gate <- function(gate, ...) UseMethod("transform_gate")

#' @export
scale_gate <- function(gate, ...) UseMethod("scale_gate")

#' @export
rotate_gate <- function(gate, ...) UseMethod("rotate_gate")

#' @export
shift_gate <- function(gate, ...) UseMethod("shift_gate")
