## ==========================================================================
## Methods for objects of type 'curv2Filter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================
## FIXME: What do summarizeFilter methods do and why are they not
## exported in the name space?
setMethod("summarizeFilter",signature("filterResult","curv2Filter"),
          function(result,filter) {
	ret = callNextMethod()
	ret$polygons = attr(result@subSet, "polygon")
	ret
})




## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("curv2Filter"),function(object) {
	msg = paste("2D curvature filter '",object@filterId,"' in dimensions ",
        paste(object@parameters, collapse=" and "), "\nwith settings:",
        "\n  bwFac=", object@bwFac, "\n  gridsize=",
        paste(object@gridsize, collapse=",", sep=""), sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})





