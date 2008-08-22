## ==========================================================================
## Methods for objects of type 'quadGate'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature(object="quadGate"),function(object)
      {
          cat("Quadrant gate '", identifier(object),
              "' with dimensions:\n", sep="")
          for(i in seq(along=object@parameters)) {
              cat("  ")
              cat(object@parameters[i])
              cat(": ")
              cat(object@boundary[i])
              cat("\n")
          }
      })



## FIXME: What is the concept of summarizeFilter methods? They are
## not exported in the namespace.
setMethod("summarizeFilter",signature("filterResult","quadGate"),
          function(result,filter) {
	ret <- callNextMethod()
	ret$populations <- levels(result@subSet)
	ret
})
