## ==========================================================================
## Methods for objects of type 'curv1Filter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================
## FIXME: What do summarizeFilter methods do and why are they not
## exported in the name space?
setMethod("summarizeFilter",signature("filterResult","curv1Filter"),
          function(result,filter) {
	ret = callNextMethod()
	ret$boundaries = attr(result@subSet, "boundaries")
	ret
})


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("curv1Filter"),function(object) {
	msg = paste("1D curvature filter '",object@filterId,"' in dimension ",
        object@parameters, "\nwith settings:",
        "\n  bwFac=", object@bwFac, "\n  gridsize=",
        paste(object@gridsize, collapse=",", sep=""), sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})








