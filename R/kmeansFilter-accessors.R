## ==========================================================================
## Methods for objects of type 'kmeansFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================



## ==========================================================================
## length method
## ---------------------------------------------------------------------------
setMethod("length",
          signature("kmeansFilter"),
          function(x) length(x@populations))

