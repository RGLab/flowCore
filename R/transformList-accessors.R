## ==========================================================================
## transformList allow us to programmatically access transformations.
## They normaly contain items of class transformMap.
## ==========================================================================






## ==========================================================================
## colnames method: This gives us the parameter names we want to transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("colnames",
          signature=signature(x="transformList"),
          definition=function(x, do.NULL=TRUE, prefix="col")
      {
          unique(sapply(x@transforms, slot, "input"))
      })



## ==========================================================================
## Concatenate two transformLists (or one list and a transformList)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("c",
          signature=signature(x="transformList"),
          definition=function(x, ..., recursive=FALSE)
      {
          ## Try to coerce everyone to a transformList first
          all.t <- lapply(list(...), as, "transformList")
          params <- c(sapply(x@transforms, slot, "output"),
                      unlist(lapply(all.t, function(x)
                                    sapply(x@transforms, slot, "output"))))
          if(length(params) > length(unique(params)))
              stop("All output parameters must be unique when combining ",
                   "transforms.", call.=FALSE)
          new("transformList", transforms=c(x@transforms,
                               unlist(sapply(all.t, slot, "transforms"))))
      })


