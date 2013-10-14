## ==========================================================================
## filterResultLists allow us to dispatch on list of filterResults as
## produced by applying a filter to a whole flowSet.
## ==========================================================================






## ==========================================================================
## Subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## to filterResultList
setMethod("[",
          signature=signature(x="filterResultList"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          if(missing(i) && missing(j)) 
              return(x)   
          if(drop)
              warning("Argument 'drop' is ignored", call.=FALSE)
          if(is.numeric(i) || is.logical(i)) {
              copy <- names(x)[i]
          } else {
              copy <- i
              i <- match(i,names(x))
          }
          if(any(is.na(copy)))
              stop("Subset out of bounds", call.=FALSE)
          if(!missing(j)){
              if(length(unique(listLen(x@.Data[i]))) !=1)
                  stop("Unequal number of populations in the",
                       "individual filterResults.\nSubsetting",
                       " not possible")
              x@.Data <- lapply(x@.Data[i], function(y) y[[j]])
          }else{
              x@.Data <- x@.Data[i]
          }
          x@frameId <- x@frameId[i]
          if(length(x@filterId)>1)
              x@filterId <- x@filterId[i]
          return(x)
      })


## to filterResult
setMethod("[[",
          signature=signature("filterResultList"),
          definition=function(x, i, j, ...)
      {
          if(length(i) != 1)
              stop("subscript out of bounds (index must have length 1)")
          if(!missing(j))
              warning("Ignoring invalid dimension", call.=FALSE)
          if(is.character(i))
              i <- match(i, x@frameId)
           if(is.na(i))
              stop("Subset out of bounds", call.=FALSE)
          return(x@.Data[[i]])
      })



## Return a more machine-readable output in form of a data.frame
setMethod("toTable",
          signature=signature(x="filterSummaryList"),
          definition=function(x, ...) {
              res <- data.frame()
              for(i in seq_along(x))
                  res <- rbind(res, data.frame(sample=names(x)[i],
                                               population=x[[i]]@name,
                                               percent=x[[i]]$p*100,
                                               count=x[[i]]@count,
                                               true=x[[i]]@true,
                                               false=x[[i]]$false,
                                               p=x[[i]]$p, q=x[[i]]$q))
              rownames(res) <- NULL
              return(res)
          })

