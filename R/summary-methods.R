## ==========================================================================
## summary methods provide useful summary statistics about an object. In
## flowCore, this will be mainly summaries about filtering opration, which
## are represented by their own class
## ==========================================================================






## ==========================================================================
## For a flowFrame, we provide summary statistics of the raw data for each
## parameter.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="flowFrame"),
          definition=function(object, ...)
      {
          apply(exprs(object), 2, summary)
      })



## ==========================================================================
## For a flowSet, we iterate over each frame and provide summaries for those.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="flowSet"),
          definition=function(object, ...)
      {
          fsApply(object, function(x) apply(exprs(x), 2, summary),
                  simplify=FALSE)
      })

## ==========================================================================
## For filters we analyse the filterResult and create a new object of class
## filterSummary.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="filter"),
          definition=function(object, result, ...)
      {
          if(missing(result))
              stop("Only resolved filters may be summarized.\nTry something ",
                   "like 'filter(myFlowFrame, myFilter)'\nbefore calling ",
                   "'summary'.", call.=FALSE)
          l <- as(result, "logical")
          true <- sum(l)
          count <- length(l)
          new("filterSummary", name=identifier(object), true=true,
              count=count, p=true/count)
      })



## ==========================================================================
## For filterReferences, we first resolve the reference to a concrete filter
## and call the next applicable method.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="filterReference"),
          definition=function(object, ...)
          summary(as(object, "concreteFilter"), ...))



## ==========================================================================
## summarize a filterResult object. There are more specific summary methods
## for the other filterResult subclasses, and this is only the fallback
## option
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="filterResult"),
          definition=function(object, ...)
      {
          summary(filterDetails(object, object@filterId)$filter,
                  object, ...)
      })



## ==========================================================================
## Summaries for logicalFilterResult objects create filterSummary objects
## with scalar values.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="logicalFilterResult"),
          definition=function(object, ...)
      {
                    filter <- filterDetails(object, identifier(object))$filter
                    identifier(filter) <- paste(identifier(filter), "+", sep="")
                    summary(filter, object)
      })



## ==========================================================================
## Summaries for multipleFilterResult objects create filterSummary objects
## with vectors
## ---------------------------------------------------------------------------
#' @export
setMethod("summary",
          signature=signature(object="multipleFilterResult"),
          definition=function(object, ...)
      {
          true <- summary(object@subSet[!is.na(object@subSet)])
          ## count <-  sum(!is.na(object@subSet))
          count <- length(object@subSet)
          name=levels(object@subSet)
          new("filterSummary", name=levels(object@subSet), true=true,
              count=count ,p=true/count)
      })



## ==========================================================================
## Summaries for manyFilterResult objects create lists of filterSummary
## objects.
## ---------------------------------------------------------------------------
#' @export
setMethod("summary",
          signature=signature(object="manyFilterResult"),
          definition=function(object, ...)
      {
          n <- names(object)
          lapply(n,function(i) summary(object[[i]]))
      })



## ==========================================================================
## summary method for filterResultLists. We print all information on the
## screen and return an object of class filterSummaryList.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="filterResultList"),
          definition=function(object, verbose=TRUE, ...)
      {
                    sums <- lapply(object, function(x){
              if(verbose)
                  cat("filter summary for frame '", x@frameId, "'\n", sep="")
              tmp <- print(summary(x), indent=1, verbose=verbose)
              if(verbose)
                  cat("\n")
              tmp
          })
          return(invisible(new("filterSummaryList", .Data=sums)))
      })



## ==========================================================================
## For subsetFilters we evaluate each component filter separately and then
## combine the results.
## --------------------------------------------------------------------------
#' @export
setMethod("summary",
          signature=signature(object="subsetFilter"),
          definition=function(object,result,...)
      {
          if(missing(result)) {
              e1 <- as(object@filters[[1]],"logical")
              e2 <- as(object@filters[[2]],"logical")
              true <- sum(e1&e2)
              ## count <- sum(e2)
              count <- length(e2)
              new("filterSummary", name=identifier(object), true=true,
                  count=count, p=true/count)
          } else {
              id <- gsub("\\+$", "", identifier(object))
              true <- sum(as(result, "logical"))
              ## count <- filterDetails(result, id)$subsetCount
              count <- length(as(result, "logical"))
              new("filterSummary", name=identifier(object), true=true,
                  count=count, p=true/count)
          }
      })



## ==========================================================================
## We want to can get parameters back OUT of a transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("summary",
          signature=signature(object="transform"),
          definition=function(object, ...)
      {
          e <- environment(object)
          x <- ls(envir=e)
          structure(lapply(x, "get", env=e), names=x)
      })
