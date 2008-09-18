



## ==========================================================================
## length method, how many populations do we have?
## ---------------------------------------------------------------------------
setMethod("length","multipleFilterResult",function(x) nlevels(x@subSet))



## ==========================================================================
## names method, the levels of the subSet factor
## ---------------------------------------------------------------------------
setMethod("names", "multipleFilterResult",function(x) levels(x@subSet))
setReplaceMethod("names", "multipleFilterResult", function(x, value)
             {
                 if(length(value) != length(levels(x@subSet)))
                     stop("Length of replacement vector doesn't match.")
                 levels(x@subSet) <- value
                 x@filterDetails[[1]]$populations <- value
                 return(x)
             })

## ==========================================================================
## subsetting method: we create a new logicalFilterResult
## ---------------------------------------------------------------------------
setMethod("[[",
          signature=signature("multipleFilterResult"),
          definition=function(x,i,j,drop=FALSE)
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


setMethod("[",
          signature=signature("multipleFilterResult"),
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


