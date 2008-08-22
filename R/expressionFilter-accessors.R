## ==========================================================================
## Methods for objects of type 'expressionFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


          
## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",
          signature("expressionFilter"),
          function(object)
      {
	msg <- paste("expression filter '", identifier(object),
                     "' evaluating the expression:\n",
                     paste(object@deparse, collapse="\n"), sep="")
	cat(msg)
	cat("\n")
	invisible(msg)
})


