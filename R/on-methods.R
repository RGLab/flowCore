## ==========================================================================
## %on% methods are basically constructors for transformFilter or
## transformList objects. They allow for an inline-syntax for transformation
## operations.
## ==========================================================================






## ===========================================================================
## Constructor for a transformFilter for a transform
## ---------------------------------------------------------------------------
setMethod("%on%",
          signature=signature(e1="filter",
                              e2="transform"),
          definition=function(e1, e2)
      {
          new("transformFilter",
              filterId=paste(identifier(e1), "on transformed values"),
              transforms=parameterTransform(e2, parameters(e1)),
              filter=e1, parameters=parameters(e1))
      })



## ===========================================================================
## Wrapper for a flowFrame. We transform on the actual raw data
## ---------------------------------------------------------------------------
setMethod("%on%",
          signature=signature(e1="transform",
                              e2="flowFrame"),
          definition=function(e1, e2)
      {
          exprs(e2) <- e1(exprs(e2))
          e2
      })



## ===========================================================================
## Constructor for a transformFilter from a parameterTransform
## ---------------------------------------------------------------------------
setMethod("%on%",
          signature=signature(e1="filter",
                              e2="parameterTransform"),
          definition=function(e1, e2)
      {
          new("transformFilter",
              filterId=paste(e1@filterId," on transformed values of ",
              paste(e2@parameters, sep=","),collapse=" "),
              transforms=e2, filter=e1,
              parameters=unique(c(e1@parameters, e2@parameters)))
})



## ===========================================================================
## Wrapper for a flowFrame. We transform on the actual raw data
## ---------------------------------------------------------------------------
setMethod("%on%",
          signature=signature(e1="parameterTransform",
                              e2="flowFrame"),
          definition=function(e1, e2)
      {
          params <- if(length(e1@parameters)) colnames(e2) else e1@parameters
          exprs(e2)[,params] <- e1(exprs(e2))[,params]
          e2
      })



## ===========================================================================
## Constructor for a transformFilter from a transformList
## ---------------------------------------------------------------------------
setMethod("%on%",
          signature=signature(e1="filter",
                              e2="transformList"),
          definition=function(e1,e2)
      {
          new("transformFilter",
              filterId=paste(e1@filterId, "on transformed values of",
              paste(colnames(e2), collapse=",")),
              transforms=e2,
              filter=e1)
      })



## ==========================================================================
## This is the real workhorse for doing %on% type transformations.
## Eventually, all %on% methods are going to call this.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%on%",
          signature=signature(e1="transformList",
                              e2="flowFrame"),
          definition=function(e1, e2)
      {
          x <- exprs(e2)
          pars <- parameters(e2)
          ranges <- range(e2)
          cN <- colnames(x)
          for(y in e1@transforms){
              if( !(y@output %in% cN) )
                  stop(y@output, "is not a variable in the flowFrame")
              x[,y@output] = y@f(x[, y@input])
              mt <- match(y@output, pars$name)
              pData(pars)[mt,c("minRange", "maxRange")] <-
                  y@f(ranges[,y@input])   
          }
          exprs(e2) <- x
          parameters(e2) <- pars
          e2
      })



## ==========================================================================
## General %on% implementation for a flowSet. Basically a wrapper around
## fsApply
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("%on%",
          signature=signature(e1="ANY",
                              e2="flowSet"),
          definition=function(e1, e2) fsApply(e2, "%on%", e1=e1))
