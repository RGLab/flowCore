## ==========================================================================
## Methods for objects of type 'complementFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("complementFilter"),
          function(object)
      {
          cat("filter '", identifier(object),
              "', the complement of\n", sep="")
          print(object@filters[[1]])
    
      })
