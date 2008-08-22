## ==========================================================================
## Methods for objects of type 'transformFilter'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## ---------------------------------------------------------------------------
setMethod("show",signature("transformFilter"),
          function(object)
      {
          cat("transformed filter '", identifier(object), "'\n", sep="")
      })



## ===========================================================================
## Constructors
## ---------------------------------------------------------------------------
setMethod("%on%",
          signature("filter", "transformList"),
          function(e1,e2)
      {
          new("transformFilter",
              filterId=paste(e1@filterId,"on transformed values of",
              paste(colnames(e2),collapse=",")),
              transforms=e2,
              filter=e1)
      })

