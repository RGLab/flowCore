## ==========================================================================
## Methods for objects of type 'norm2Filter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
##  summarize the results of a norm2Filter operation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## FIXME: What do summarizeFilter methods do and why are they not
## exported in the name space?
setMethod("summarizeFilter",signature("filterResult","norm2Filter"),
          function(result,filter) {
              ret = callNextMethod()
	ret$cov    = attr(result@subSet,'cov')
	ret$center = attr(result@subSet,'center')
	ret$radius = attr(result@subSet,'radius')
	ret
})


## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
          signature(object="norm2Filter"),
          function(object)
      {
          cat(ifelse(length(object@transformation), "transformed", ""),
              "norm2Filter '", identifier(object),
              "' in dimensions ", sep="")
          cat(paste(object@parameters, sep="", collapse=" and "),
              "with parameters:\n")
          cat("  method:", object@method, "\n")
          cat("  scale.factor:", object@scale.factor, "\n")
          cat("  n:", object@n, "\n")
          cat("\n")
      })
