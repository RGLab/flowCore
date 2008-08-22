## ==========================================================================
## Methods for objects of type 'sampleFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("sampleFilter"),function(object) {
	msg = paste("sample filter '", object@filterId,
        "' returning objects with ", object@size," rows", sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})
