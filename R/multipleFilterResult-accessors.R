## ==========================================================================
## multipleFilterResults create a multiple non-overlapping populations
## ==========================================================================


## ==========================================================================
## Subsetting to the individual populations. We create a new
## logicalFilterResult
## ---------------------------------------------------------------------------
#' @export
setMethod("[[",
          signature=signature(x="multipleFilterResult"),
          definition=function(x, i, j, drop=FALSE)
      {
          if(length(i)!=1)
              stop("subscript out of bounds (index must have length 1)")
          if(is.numeric(i))
              i <- names(x)[i]
          if(is.na(i) || !i %in% names(x))
              stop("filter index out of bounds", call.=FALSE)
          filterDetails <- x@filterDetails
          filterDetails$population <- i
          new("logicalFilterResult",subSet=(x@subSet==i & !is.na(x@subSet)),
              filterDetails=filterDetails,
              frameId=x@frameId, filterId=paste(x@filterId,"[[",i,"]]",
                                 sep=""))
      })



## ==========================================================================
## Subsetting to the a (smaller) multipleFilterResult 
## ---------------------------------------------------------------------------
#' @export
setMethod("[",
          signature=signature(x="multipleFilterResult"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          if(is.numeric(i))
              i <- names(x)[i]
          if(any(is.na(i)) || !all(i %in% names(x)))
              stop("filter index out of bounds", call.=FALSE)
          filterDetails <- x@filterDetails
          filterDetails$population <- i
          x@subSet <- factor(x@subSet[x@subSet %in% i & !is.na(x@subSet)])
          x@filterDetails <- filterDetails
          return(x)
      })


