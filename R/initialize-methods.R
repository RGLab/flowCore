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

setMethod("initialize", "flowFrame",
          function(.Object,
                   exprs,
                   parameters,
                   description=list(note="empty"), ...){
              if (missing(parameters))
                parameters <- flowDefaultADF(exprs)
              
              callNextMethod(.Object,
                             exprs=exprs,
                             parameters=parameters,
                             description=description,
                             ...)
          })

