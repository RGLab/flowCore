## ==========================================================================
## Methods for objects of type 'polygoneGate'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
          signature(object="polygonGate"),
          function(object)
      {
          nb <-  nrow(object@boundaries)
          cat("Polygonal gate '", identifier(object) ,"' with ",
              ifelse(all(is.na(object@boundaries)), 0, nb),
              " vertices in dimensions ", sep="")
          cat(paste(object@parameters, sep="", collapse=" and "))
          cat("\n")
      })
