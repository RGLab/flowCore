## ==========================================================================
## Show the length of a (usually list-like) object. For many flowCore
## gating objects this returns the number of sub-populations.
## ==========================================================================






## ==========================================================================
## Number of populations (i.e., k) defined by a kmeansFilter
## ---------------------------------------------------------------------------
#' @export
setMethod("length",
          signature=signature(x="kmeansFilter"),
          definition=function(x) length(x@populations))



## ==========================================================================
## Essentially the number of populations (in a multipleFilterResult)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("length",
          signature=signature(x="filterSummary"),
          definition=function(x) length(x@name))



## ==========================================================================
## The number of frames in a flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("length",
          signature=signature(x="flowSet"),
          definition=function(x) length(sampleNames(x)))



## ==========================================================================
## For logicalFilterResults we only have a single population, hence the
## length is always 1
## ---------------------------------------------------------------------------
#' @export
setMethod("length",
          signature=signature(x="logicalFilterResult"),
          definition=function(x) 1)



## ==========================================================================
## For manyFilterResults, length corresponds to the number of
## sub-populations defined by the various component filters
## ---------------------------------------------------------------------------
#' @export
setMethod("length",
          signature=signature(x="manyFilterResult"),
          definition=function(x) ncol(x@subSet))



## ==========================================================================
## For multipleFilterResults, the number of sub-populations
## ---------------------------------------------------------------------------
#' @export
setMethod("length",
          signature=signature(x="multipleFilterResult"),
          definition=function(x) nlevels(x@subSet))
