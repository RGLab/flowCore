## ==========================================================================
## Methods for objects of type 'intersectFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("intersectFilter"),
          function(object)
      {
          cat("filter '", identifier(object),
              "'\nthe intersection between the ", length(object@filters),
              " filters\n\n", sep="")
          for(i in 1:length(object@filters)){
              print(object@filters[[i]])
              cat("\n")
          }
      })



