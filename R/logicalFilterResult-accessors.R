## ==========================================================================
## logicalFilterResults create a single population (and the complement of it)
## ==========================================================================






## ==========================================================================
## Subset a logicalFilterResult by filterId. We treat i=1 as the population
## in the filter and 1=2 as the complement of that. Values >2 are not allowed
## --------------------------------------------------------------------------
setMethod("[[",
          signature=signature(x="logicalFilterResult"),
          definition=function(x, i, j, drop=FALSE)
      {
          if(drop)
              warning("Argument 'drop' is ignored", call.=FALSE)
          if(!missing(j))
              warning("Ignoring invalid dimension", call.=FALSE)
          if(is.numeric(i) || is.logical(i)) {
              copy <- names(x)[i]
          }else {
              copy <- i
              i <- match(i,names(x))
          }
          if(is.na(copy))
              stop("Subset out of bounds", call.=FALSE)
          if(i==1){
              x
          }else{
              tmp <- x
              tmp@subSet <- !x@subSet
              fd <- filterDetails(x)
              id <- identifier(fd$filt)
              fd$filter <- !fd$filt
              fd$populations <- rev(fd$populations)
              filterDetails(tmp, id) <- fd
              tmp
          }   
      })
