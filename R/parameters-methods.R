## ==========================================================================
## Parameters in the context of flow cytometry are the measurement
## parameters of the experiment, i.e., the columns in the flowFrame.
## These methods return or set parameters for the various objects in
## flowCore.
## ==========================================================================






## ==========================================================================
## Default for everything that inherits from filter. Since no parameters
## are defined at that level we return an empty charactor vector
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="filter"),
          definition=function(object) character(0))



## ==========================================================================
## For paramterFilters we can extract things directly from the parameters
## slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="parameterFilter"),
          definition=function(object){
            tmp <- unlist(sapply(object@parameters, as, "character"))
            if(all(sapply(tmp, is.character)))
              tmp
            else
              NULL
          })



setReplaceMethod("parameters", 
                 signature=signature(object="parameterFilter",
                   value="character"),
                 definition=function(object, value)
                 {
                   value <- sapply(value, unitytransform)   
                   object@parameters <- new("parameters",
                                            .Data=value)
                   object
                 })


setReplaceMethod("parameters", 
                 signature=signature(object="parameterFilter",
                   value="transform"),
                 definition=function(object, value)
                 {
                   object@parameters <- new("parameters",
                                            .Data=value)
                   object
                 })


setReplaceMethod("parameters", 
                 signature=signature(object="parameterFilter",
                   value="list"),
                 definition=function(object, value)
                 {
                   value <- unlist(value)
                   chars <- sapply(value, is.character)
                   if(!is.null(chars) && any(chars))
                     value[chars] <- sapply(value[chars],
                                              unitytransform)
                   if(!all(sapply(value, is, "transform")))
                     stop("Items in the list need to be characters or ",
                          "trasnforms.", call.=FALSE)
                   object@parameters <- new("parameters",
                                            .Data=value)
                   object
                 })


## ==========================================================================
## For transformations
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="transform"),
          definition=function(object){
            tmp <- unlist(sapply(object@parameters, as, "character"))
            tmp
          })

setMethod("parameters",
          signature=signature(object="ratio"),
          definition=function(object) NA)


setMethod("parameters",
          signature=signature(object="nullParameter"),
          definition=function(object) NA) 


setReplaceMethod("parameters", 
                 signature=signature(object="singleParameterTransform",
                   value="transform"),
                 definition=function(object, value)
                 {
                   object@parameters <- value
                   object
                 })

setReplaceMethod("parameters", 
                 signature=signature(object="dg1polynomial",
                   value="transform"),
                 definition=function(object, value)
                 {
                   object@parameters <- new("parameters",
                                            .Data=list(value))
                   object
                 })


setReplaceMethod("parameters", 
                 signature=signature(object="singleParameterTransform",
                   value="character"),
                 definition=function(object, value)
                 {
                   value <- unitytransform(value)   
                   object@parameters <- value
                   object
                 })


## ==========================================================================
## For setOperationFilters we return the parameters of every subfilter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="setOperationFilter"),
          definition=function(object)
          unique(unlist(lapply(object@filters, parameters))))



## ==========================================================================
## For a filterReference we first resolve the reference and use the next
## available method for the concrete filter object.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature(object="filterReference"),
          function(object) parameters(as(object, "concreteFilter")))



## ==========================================================================
## For filterResults we can get the parameters from the filter object in
## the fitlerDetails slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="filterResult"),
          definition=function(object) filterDetails(object)$parameters)



## ==========================================================================
## For manyFilterResults we need to parse all items in the filterDetails
## list
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="manyFilterResult"),
          definition=function(object)
          unique(sapply(filterDetails(object),function(x) x$parameters)))



## ==========================================================================
## In a flowFrame we can get things directly from the parameters slot. This
## returns the complete annotatedDataFrame unless names=TRUE
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="flowFrame"),
          definition=function(object, names=FALSE)
      {
          if(!names)
              object@parameters
          else
              as.character(parameters(object)[["name"]])
      })



## ==========================================================================
## For filterResultsList we have to parse through the individual
## filterResults
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="filterResultList"),
          definition=function(object) lapply(object, parameters))



## ==========================================================================
## Replacement of parameters in a flowFrame is only possible with
## another annotatedDataFrame.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("parameters",
                 signature=signature(object="flowFrame",
                 value="AnnotatedDataFrame"),
                 definition=function(object, value){
                     if(!all(c("name", "desc", "range", "minRange",
                               "maxRange") %in% varLabels(value)))
                         stop("varLabels of this AnnotatedDataFrame don't ",
                              "match the specifications", call.=FALSE)
                     if(!all(colnames(exprs(object)) ==  value$name))
                         stop("parameter names don't match colnames of the ",
                              "exprs matrix", call.=FALSE)
                     object@parameters <- value
                     return(object)
                 })



## ==========================================================================
## In a parameterTransform we can get things directly from the parameters
## slot.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("parameters",
          signature=signature(object="parameterTransform"),
          definition=function(object) object@parameters)

