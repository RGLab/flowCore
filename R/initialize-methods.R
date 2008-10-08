## ===========================================================================
## If no usuful parameters slot is provided for a flowFrame we create one
## ---------------------------------------------------------------------------
flowDefaultADF <- function(exprs) {
    ## contract: 'object' is matrix-like, with dim, rownames, colnames
    ## methods. Returns AnnotatedDataFrame with appropriate dimensions.
    vm <- data.frame(labelDescription=c(name="Name of Parameter",
                       desc="Description of Parameter",
                       range="Range of Parameter",
                       minRange="Minimum Parameter Value after Transformation",
                       maxRange="Maximum Parameter Value after Transformation"))
    pd <- data.frame(name=colnames(exprs), desc=colnames(exprs),
                     range=apply(exprs, 2, max, na.rm=TRUE),
                     minRange=apply(exprs, 2, min, na.rm=TRUE),
                     maxRange=apply(exprs, 2, max, na.rm=TRUE))
    new("AnnotatedDataFrame", pd, vm)
}

setMethod("initialize",
          signature=signature(.Object="flowFrame"),
          definition=function(.Object, exprs, parameters,
                              description=list(note="empty"), ...)
      {
          if (missing(parameters))
              parameters <- flowDefaultADF(exprs)
          
          callNextMethod()
      })



## ===========================================================================
## Make sure that the old filter class constructors still work with
## the new 'parameter' slot. We archive this by populating the parameters
## slot using the assignment method which can deal with different inputs
## ---------------------------------------------------------------------------
setMethod("initialize",
          signature=signature(.Object="parameterFilter"),
          definition=function(.Object, parameters, ...)
          {
            parameters(.Object) <- parameters
            callNextMethod()
          })



## ===========================================================================
## For transforms we need to initialize both the parameter slot (if there is
## one) and also the function in the .Data slot.
## ---------------------------------------------------------------------------
## Parameters for all single parameter transformations
setMethod("initialize",
          signature=signature(.Object="singleParameterTransform"),
          definition=function(.Object, parameters, ...)
          {
            parameters(.Object) <- parameters
            callNextMethod()
          })

## this is a multi parameter transform, we need to set up the parameters
## explicitely
setMethod("initialize",
          signature=signature(.Object="dg1polynomial"),
          definition=function(.Object, parameters, ...)
          { 
            parameters(.Object) <- parameters
            .Object@.Data <- function(...)
            {   
                temp <- sapply(list(...), cbind)
                coeff <- .Object@a
                res <- 0
                for(i in seq_len(ncol(temp)))
                    res <- res + temp[,i,drop=FALSE]*coeff[i]
                res <- res+.Object@b
            }
            callNextMethod()
          })

## this is a multi parameter transform, we need to set up the parameters
## explicitely
setMethod("initialize",
          signature=signature(.Object="ratio"),
          definition=function(.Object, ...)
      { 
          .Object@.Data <- function(x, y) x/y
          callNextMethod()
      })


