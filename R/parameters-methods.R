## ==========================================================================
## Parameters in the context of flow cytometry are the measurement
## parameters of the experiment, i.e., the columns in the flowFrame.
## These methods return or set parameters for the various objects in
## flowCore.
## ==========================================================================





#' Obtain information about parameters for flow cytometry objects.
#' 
#' 
#' Many different objects in \code{flowCore} are associated with one or more
#' parameters. This includes \code{\linkS4class{filter}},
#' \code{\linkS4class{flowFrame}} and \code{\linkS4class{parameterFilter}}
#' objects that all either describe or use parameters.
#' 
#' 
#' @name parameters-methods
#' @aliases parameters parameters,flowFrame,missing-method
#' parameters,filter-method parameters,filterResult-method
#' parameters,parameterFilter-method parameters,setOperationFilter-method
#' parameters,filterReference-method parameters,parameterTransform-method
#' parameters,flowFrame-method parameters<-
#' parameters<-,flowFrame,AnnotatedDataFrame-method
#' parameters<-,dg1polynomial,transform-method
#' parameters<-,parameterFilter,character-method
#' parameters<-,parameterFilter,list-method
#' parameters<-,parameterFilter,transform-method
#' parameters<-,singleParameterTransform,character-method
#' parameters<-,singleParameterTransform,transform-method
#' parameters,nullParameter-method parameters,ratio-method
#' parameters,transform-method
#' @docType methods
#' @usage parameters(object, \dots)
#' 
#' @param object Object of class \code{\linkS4class{filter}},
#' \code{\linkS4class{flowFrame}} or \code{\linkS4class{parameterFilter}}.
#' @param \dots Further arguments that get passed on to the methods.
#' @return
#' 
#' When applied to a \code{flowFrame} object, the result is an
#' \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#' describing the parameters recorded by the cytometer. For other objects it
#' will usually return a vector of names used by the object for its
#' calculations.
#' 
#' @section Methods:
#' \describe{
#' 
#' \item{parameters(object = "filter")}{Returns for all objects that inherit from
#' \code{\linkS4class{filter}} a vector of parameters on which a gate is
#' defined.}
#' 
#' \item{parameters(object = "parameterFilter")}{see above}
#' 
#' \item{parameters(object = "setOperationFilter")}{see above}
#' 
#' \item{parameters(object = "filterReference")}{see above}
#' 
#' \item{parameters(object = "flowFrame")}{Returns an
#' \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
#' containing detailed descriptions about the measurement parameters of the
#' \code{\linkS4class{flowFrame}}. For \code{\linkS4class{flowFrame}} objects
#' there also exists a replacement method.}
#' 
#' }
#' @author B. Ellis, N. Le Meur, F. Hahne
#' @keywords methods
#' @examples
#' 
#'  samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
#'  parameters(samp)
#'  print(samp@parameters@data)
#' 

## ==========================================================================
## Default for everything that inherits from filter. Since no parameters
## are defined at that level we return an empty charactor vector
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="filter"),
          definition=function(object) character(0))



## ==========================================================================
## For paramterFilters we can extract things directly from the parameters
## slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="parameterFilter"),
          definition=function(object){
            tmp <- unlist(sapply(object@parameters, as, "character"))
            if(all(sapply(tmp, is.character)))
              tmp
            else
              NA
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
                                            .Data=list(value))
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
#' @export
setMethod("parameters",
          signature=signature(object="transform"),
          definition=function(object){
            tmp <- unlist(sapply(object@parameters, as, "character"))
            tmp
          })

#' @export
setMethod("parameters",
          signature=signature(object="ratio"),
          definition=function(object) NA)

setMethod(
    "parameters",
    signature = signature(object = "ratiotGml2"),
    definition = function(object) NA
)
 
  
# J. Spidlen: Oct 22, 2013
# The unlist(sapply(object@parameters, as, "character")) that was normally used for 
# singleParameterTransform gives:
# Error in as.vector(x, "list") : cannot coerce type 'closure' to vector of type 'list'
# Returning object@transformationId instead seems to work just fine. 
#' @export
setMethod("parameters",
          signature=signature(object="singleParameterTransform"),
          definition=function(object) {
              object@transformationId
		  })

#' @export
setMethod("parameters",
          signature=signature(object="nullParameter"),
          definition=function(object) NA)


#' @export
setMethod("parameters",
          signature=signature(object="transformReference"),
          definition=function(object){
          while(is(object, "transformReference"))
            object <- eval(object)
          parameters(object)
        })


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
                 signature=signature(object="dg1polynomial",
                   value="transform"),
                 definition=function(object, value)
                 {
                   object@parameters <- new("parameters",
                                            .Data=list(value))
                   object
                 })

setReplaceMethod("parameters", 
                 signature=signature(object="dg1polynomial",
                   value="parameters"),
                 definition=function(object, value)
                 {
                   object@parameters <- value
                   object
                 })

setReplaceMethod("parameters", 
                 signature=signature(object="dg1polynomial",
                   value="character"),
                 definition=function(object, value)
                 {
                      value <- new("parameters", .Data=list(unitytransform(value))) 
                      object@parameters <- value 
                                            
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
## For compensations
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="compensation"),
          definition=function(object){
            tmp <- unlist(sapply(object@parameters, as, "character"))
            tmp
          })



## ==========================================================================
## For setOperationFilters we return the parameters of every subfilter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="setOperationFilter"),
          definition=function(object)
          unique(unlist(lapply(object@filters, parameters))))



## ==========================================================================
## For a filterReference we first resolve the reference and use the next
## available method for the concrete filter object.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature(object="filterReference"),
          function(object) parameters(as(object, "concreteFilter")))



## ==========================================================================
## For filterResults we can get the parameters from the filter object in
## the fitlerDetails slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="filterResult"),
          definition=function(object)
          unique(unlist(sapply(filterDetails(object),
                               function(x) if(is.list(x)) x$parameters))))



## ==========================================================================
## For manyFilterResults we need to parse all items in the filterDetails
## list
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="manyFilterResult"),
          definition=function(object)
          unique(unlist(sapply(filterDetails(object),
                               function(x) x$parameters))))



## ==========================================================================
## In a flowFrame we can get things directly from the parameters slot. This
## returns the complete annotatedDataFrame unless names=TRUE
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
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
#' @export
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
#' @export
setMethod("parameters",
          signature=signature(object="parameterTransform"),
          definition=function(object) object@parameters)



# ==========================================================================
## In a normalization we can get things directly from the parameters
## slot.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("parameters",
          signature=signature(object="normalization"),
          definition=function(object) object@parameters)
