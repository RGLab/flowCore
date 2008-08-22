## ==========================================================================
## Methods for objects of type 'kmeansFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================
## FIXME: What do summarizeFilter methods do and why are they not
## exported in the name space?
setMethod("summarizeFilter",
          signature("filterResult","kmeansFilter"),
          function(result,filter)
      {
          ret <- callNextMethod()
          ret$populations = filter@populations
          ret
      })


## ==========================================================================
## length method
## ---------------------------------------------------------------------------
setMethod("length",
          signature("kmeansFilter"),
          function(x) length(x@populations))


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",
          signature("kmeansFilter"),
          function(object)
      {
          msg <- paste("k-means filter '", object@filterId,
                       "' in dimension ", object@parameters[1],
                       "\nwith ", length(object), " populations (",
                       paste(object@populations, collapse=","),
                       ")", sep="")
          cat(msg)
          cat("\n")
          invisible(msg)
      })
