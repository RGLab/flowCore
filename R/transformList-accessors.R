## ==========================================================================
## transformList allow us to programmatically access transformations.
## They normaly contain items of class transformMap.
## ==========================================================================

## ==========================================================================
## colnames method: This gives us the parameter names we want to transform
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",
          signature("transformList"),
          function(x, do.NULL=TRUE, prefix="col")
      {
          unique(sapply(x@transforms, slot, "input"))
      })



## ==========================================================================
## %on% operators: this is the real workhorse for doing %on% type
## transformations. Eventually, all %on% methods are going to call this
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%on%",
          signature("transformList", "flowFrame"),
          function(e1, e2) {
              x <- exprs(e2)
              pars <- parameters(e2)
              ranges <- range(e2)
              cN <- colnames(x)
              for(y in e1@transforms){
                  if( !(y@output %in% cN) )
                      stop(y@output, "is not a variable in the flowFrame")
                  x[,y@output] = y@f(x[, y@input])
                  mt <- match(y@output, pars$name)
                  pData(pars)[mt,c("minRange", "maxRange")] <-
                      y@f(ranges[,y@input])   
              }
              exprs(e2) <- x
              parameters(e2) <- pars
              e2
          })



## ==========================================================================
## coerce from a list to a transformLis
setAs(from="list", to="transformList",
      def=function(from) new("transformList", transforms=from))



## ==========================================================================
## concatenate two transformLists (or one list and a transformList)
setMethod("c",
          signature("transformList"),
          function(x, ..., recursive=FALSE)
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



## ==========================================================================
## print details about the transformMaps in a transformList
setMethod("show",
          signature("transformMap"),
          function(object) cat("transformMap for parameter '",
                               object@input, "' mapping to '",
                               object@output, "'\n", sep=""))
