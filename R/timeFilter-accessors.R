## ==========================================================================
## Methods for objects of type 'timeFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("timeFilter"),function(object)
          {
            msg <- paste("time filter '",object@filterId,
                         "' with settings:\n  bandwidth=",
                         object@bandwidth, sep="")
            cat(msg)
            if(length(object@binSize))
              cat("\n  binSize=", object@binSize, sep="")
            if(length(object@timeParameter))
              cat("\n  timeParameter=", object@timeParameter, sep="")
            cat("\n")
            invisible(msg)
          })


