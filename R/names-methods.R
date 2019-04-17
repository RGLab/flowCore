## ==========================================================================
## Show the name of an object, or of the items contained in a (usually
## list-like) object. In flowCore, names are often used for subpopulations.
## Note that parameter names are accessed by the parameters method.
## ==========================================================================


## ==========================================================================
## The identifiers of the flowFrames (i.e., the sampleNames of the flowSet)
## for which the filterResults have been computed.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("names",
          signature=signature(x="filterResultList"),
          definition=function(x) x@frameId)



## ==========================================================================
## The names of the populations if the summary was computed for a
## multipleFilterResult, for a logicalFilterResult the name of the input
## filter.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("names",
          signature=signature(x="filterSummary"),
          definition=function(x) x@name)



## ==========================================================================
## This return a pretified version of the parameter names, including the
## parameter description if present
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("names",
          signature=signature(x="flowFrame"),
          definition=function(x)
      {
          cn <- colnames(x)
          fn <- featureNames(x)
          if(length(fn) == length(cn)) {
              cn <- paste("<", cn, ">", sep="")
              for(i in seq(along=fn)) {
                  if(!is.na(fn[i]) && fn[i]!="")
                      cn[i] <- paste(cn[i],fn[i])
              }
          }
          cn
      })



## ==========================================================================
## The names of the two populations created by the filter. For
## logicalFilterResults we always have one population in the filter
## and the complement.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("names",
          signature=signature(x="logicalFilterResult"),
          definition=function(x) paste(x@filterId, c("+", "-"), sep=""))



## ==========================================================================
## The names of the individual sub-populations (i.e., the colnames of the
## subSet matrix.
## ---------------------------------------------------------------------------
#' @export
setMethod("names",
          signature=signature(x="manyFilterResult"),
          definition=function(x) c("rest", colnames(x@subSet)))



## ==========================================================================
## The names of the individual sub-populations (i.e., the levels of the
## subSet factor). The replacement method simply changes the factor levels.
## ---------------------------------------------------------------------------
#' @export
setMethod("names",
          signature=signature(x="multipleFilterResult"),
          definition=function(x) levels(x@subSet))

setReplaceMethod("names",
                 signature=signature(x="multipleFilterResult",
                                     value="ANY"),
                 definition=function(x, value)
             {
                 if(length(value) != length(levels(x@subSet)))
                     stop("Length of replacement vector doesn't match.")
                 levels(x@subSet) <- value
                 x@filterDetails[[1]]$populations <- value
                 return(x)
             })
