## ==========================================================================
## Methods for objects of type 'unionFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("unionFilter"),
          function(object)
      {
          cat("filter '", identifier(object),
              "'\nthe union of the ", length(object@filters),
              " filters\n\n", sep="")
          for(i in 1:length(object@filters)){
              print(object@filters[[i]])
              cat("\n")
          }
      })


