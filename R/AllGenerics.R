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
setGeneric("spillover", function(x,...){
  tryCatch(standardGeneric("spillover"),
      error = function(e){
        if(grepl("unable to find an inherited method for function 'spillover' for signature.*(flowSet|ncdfFlowSet)", e$message)){
          stop("The flowSet spillover method has been moved to the flowStats package.
               Please library(flowStats) first.")
        }else{
          stop(e)
        }
      })
})
#' @export
setGeneric("compensate", 
            function(x, spillover, ...)
            standardGeneric("compensate"))
#' @export
setGeneric("decompensate",
			  function(x, spillover, ...)
			  	standardGeneric("decompensate"))

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

