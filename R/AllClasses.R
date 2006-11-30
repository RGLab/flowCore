require("Biobase")

##  Classes definition - last changes Oct 27, 2006
##  pdh: redefined RectangleGate so that it correctly uses min and max instead of
##      coordinates of the corners.

## ===========================================================================
##  ~cytoFrame/FCS no gating slot / no well Slot?
## ---------------------------------------------------------------------------
## A container for flow cytometry measurements with slots exprs, description
## and well. Exprs contains measurement values, description contains 
## information from file headers of FCS file and well contains well position
## on microtiter plate from experiment.  
## 
## ---------------------------------------------------------------------------
setClass("flowFrame",                
  representation(exprs="matrix",
                 description="character"),
  prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
                 description=c(note="empty")),
 validity=function(object){
   msg <- TRUE
   if(!is.matrix(object@exprs))
      msg <- "\nslot 'exprs' must be matrix"
   if(!is.character(object@description))
      msg <- "\nslot 'description' must be character vector"
})

## ===========================================================================
##  ~cytoSet
## ---------------------------------------------------------------------------
## A collection of several cytoFrames making up one experiment. Slots 
## frames, phenoData, colnames. Frames contains the cytoFrame objects,
## phenoData the experiment meta data and colnames the channel names. 
## ---------------------------------------------------------------------------
setClass("flowSet",                   
  representation(frames="environment",
                 phenoData="AnnotatedDataFrame",
                 colnames="character"),
  prototype=list(frames=new.env(),
                 phenoData=new("AnnotatedDataFrame",
                   data=data.frame(name=I(character(0))),
                   varMetadata=data.frame(labelDescription="Name in frame",row.names="name")),
                 colnames=character(0)),
  validity=function(object){
    nc <- length(colnames(object))
	if(!(is(object@phenoData, "AnnotatedDataFrame") && 
		is(object@colnames, "character") &&
		is(object@frames, "environment"))) {
			warning("An element of this object is of the wrong type.")
			return(FALSE)
		}
	if(!("name" %in% colnames(pData(object@phenoData)))) {
		warning("phenoData has no name column")
		return(FALSE)
	}
	if(any(is.na(match(object@phenoData$name,ls(object@frames,all.names=TRUE))))) {
		warning("Some names given in phenoData do not exist in this set.")
		return(FALSE)
	}
	
	if(!all(sapply(object@phenoData$name,function(i) {
		x = get(i,env=object@frames)
		if(ncol(exprs(x)) != nc || !is.null(colnames(x))) FALSE else TRUE
	}))) {
		warning("Some items identified in set are either have the wrong dimension or type.")
		return(FALSE)
	}
	return(TRUE)
  })



## ===========================================================================
##  Transformation
## ---------------------------------------------------------------------------
## Definition of the a transformation
## 
## 
## ---------------------------------------------------------------------------
setClass("transformation", 
         representation(transformationId="character",
                        parameters="character"), 
         prototype=prototype(transformationId="ANY",
                        parameters="character"), 
         contains="VIRTUAL",
         
         validity=function(object){
             msg <- TRUE
             if(!is.character(object@transformationId) ||
                length(object@transformationId)!=1)
               msg <- "\nslot 'name' must be character vector of length 1"
             if(!is.character(object@parameters))
             msg <- "\nslot 'parameters' must be character vector"
             test <- matrix(1:length(object@parameters), ncol=length(object@parameters))
             colnames(test) <- object@parameters
             return(msg)
         })

## ===========================================================================
##  Linear Transformation
## ---------------------------------------------------------------------------
## Definition of the a linear transformation
## ---------------------------------------------------------------------------
setClass("linearTransformation",
         contains="transformation",
         representation(a="numeric",b="ANY"), 
         prototype=prototype(transformationId="linear-Transformation",a=1), 
         validity=function(object){
             msg <- TRUE
             if(!is.numeric(object@a))
               msg <- "\nslot 'a' must be numeric" 
             if (object@b && !is.numeric(object@b))
                msg <- "\nslot 'b' must be numeric"
             return(msg)
         })
## ===========================================================================
## Quadratic Transformation
## ---------------------------------------------------------------------------
## Definition of the a quadratic transformation
## ---------------------------------------------------------------------------
setClass("quadraticTransformation",
         contains="transformation",
         representation(a="numeric",b="numeric",c="numeric"), 
         prototype=prototype(transformationId="quadratic-Transformation",a=1,b=0,c=0), 
         validity=function(object){
             msg <- TRUE
             if(!is.numeric(object@a) || !is.numeric(object@b) || !is.numeric(object@c))
               msg <- "\nslot 'a', 'b', 'c' must be numeric" ##one of them might be missing?
             return(msg)
         })
## ===========================================================================
## Natural Logarithm Transformation
## ---------------------------------------------------------------------------
## Definition of the a  Natural Logarithm transformation
## ---------------------------------------------------------------------------
setClass("lnTransformation",
         contains="transformation",
         representation(r="numeric",d="numeric"), 
         prototype=prototype(transformationId="LN-Transformation",r=1,d=1), 
         validity=function(object){
             msg <- TRUE
             if((!is.numeric(object@r) & object@r >0) || (!is.numeric(object@d)& object@d >0 ))
               msg <- "\nslot 'r', 'd' must be positive number" ##one of them might be missing?
             return(msg)
         })

## ===========================================================================
## LOG Transformation
## ---------------------------------------------------------------------------
## Definition of the a Log transformation
## ---------------------------------------------------------------------------
setClass("logTransformation",
         contains="transformation",
         representation(r="numeric",d="numeric",logbase="numeric"), 
         prototype=prototype(transformationId="Log-Transformation",r=1,d=1,logbase=10), 
         validity=function(object){
             msg <- TRUE
             if((!is.numeric(object@r) & object@r >0) || (!is.numeric(object@d)& object@d >0 ))
               msg <- "\nslot 'r', 'd' must be positive number" ##one of them might be missing?
             if(!is.numeric(object@logbase) & object@d >1 )
               msg <- "\nslot 'r', 'd' must be number greater than 1" ##one of them might be missing?
             return(msg)
         })
## ===========================================================================
## HyperLOG Transformation
## ---------------------------------------------------------------------------
## Definition of the a HyperLOG transformation
## ---------------------------------------------------------------------------
setClass("hyperlogTransformation",
         contains="transformation",
         representation(r="numeric",d="numeric",b="numeric"), 
         prototype=prototype(transformationId="HyperLog-Transformation",r=1000,d=3,b=35), 
         validity=function(object){
             msg <- TRUE
             if((!is.numeric(object@b) & object@b >0))
                msg <- "\nslot 'b' )hyperlog) must be positive number"
             if(any(!is.numeric(object@r) & object@r >0) || (!is.numeric(object@d)& object@d >0 ))
               msg <- "\nslot 'r', 'd' must be positive number" ##one of them might be missing?
             return(msg)
         })

## ===========================================================================
## Biexponential Transformation
## ---------------------------------------------------------------------------
## Definition of the general biexponential transformation
## ---------------------------------------------------------------------------
setClass("biexponentialTransformation",
         contains="transformation",
         representation(a="numeric",b="numeric",c="numeric",d="numeric",f="numeric",w="numeric",tol="numeric","maxit"="integer"), 
         prototype=prototype(transformationId="Biexponential-Transformation",a=.5,b=1,c=.5,d=1,f=0,w=0,
           tol=.Machine$double.eps^0.25,maxit=as.integer(5000)), 
         validity=function(object){
           msg <- TRUE
         })


## ===========================================================================
## Virtual filter/subsetting?
## ---------------------------------------------------------------------------
## An object describing a selection applied to a data matrix. Consist of
## a functions that return logical vectors subsetting the data
## ---------------------------------------------------------------------------
setClass("filter", 
         representation(filterId="character",
                        parameters="ANY"), 
         prototype=prototype(filterId="character",
                        parameters="character"), 
         contains="VIRTUAL",
         
         validity=function(object){
             msg <- TRUE
             if(!is.character(object@filterId) ||
                length(object@filterId)!=1)
               msg <- "\nslot 'name' must be character vector of length 1"
             if(!is.character(object@parameters))
             msg <- "\nslot 'parameters' must be  vector"
            
             test <- matrix(1:length(object@parameters),
                            ncol=length(object@parameters))
             colnames(test) <- object@parameters
             return(msg)
         })
## ===========================================================================
## filterSet
## ---------------------------------------------------------------------------
## An object describing a set of individual gating functions to subset
## data, possibly in several dimensions.
## ---------------------------------------------------------------------------


## ===========================================================================
## Rectangular gate
## ---------------------------------------------------------------------------
setClass("rectangleGate",
         contains="filter",
         representation(min="numeric",
                        max="numeric"),
         ##minimum and maximum values for the coordinates
         prototype=list(filterId="Rectangle Gate",
           min=0,max=Inf)
         )

## ===========================================================================
## Polygon gate
## ---------------------------------------------------------------------------
setClass("polygonGate",
         representation(boundaries="matrix"),
         contains="filter",
         prototype=list(filterId="ALL", boundaries=matrix(ncol=2, nrow=3)),

         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@boundaries) ||
                nrow(object@boundaries)<2)
               msg <- "\nslot 'boundaries' must be character vector longer than 2"
             return(msg)
         })


## ===========================================================================
## Polytope gate
## ---------------------------------------------------------------------------
setClass("polytopeGate",
         representation(boundaries="matrix"),
         contains="filter",
         prototype=list(filterId="ALL", boundaries=matrix(ncol=2, nrow=2)), 

         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@boundaries) ||
                nrow(object@boundaries)<2)
               msg <- "\nslot 'boundaries' must be character vector longer then 2"
             return(msg)
         })


## ===========================================================================
## Ellipsoid gate
## ---------------------------------------------------------------------------
## 
## 
## ---------------------------------------------------------------------------
setClass("ellipsoidGate",
         representation(focus="matrix",
                        distance="numeric"),
         contains="filter",
         prototype=list(filterId="ALL", focus=matrix(ncol=2, nrow=2)),

         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@focus) ||
                nrow(object@focus)< 2)
               msg <- "\nslot 'focus' must be matrix of at least 2 rows"
             if(!is.numeric(object@distance) ||
                length(object@distance)!=1)
               msg <- "\nslot 'distance' must be numeric vector of length 1"
             return(msg)
         })


## ===========================================================================
## norm2Filter (adapted from prada)
## ---------------------------------------------------------------------------
## 
## 
## ---------------------------------------------------------------------------
setClass("norm2Filter",
         ## the slot method holds the method argument to fitNorm2
         ## the slot scale.factor holds the scalefac argument to fitNorm2
         ## transformation holds a list of length giving transformations, if applicable that are
         ## applied to the data before gating
         representation(method="character",
                        scale.factor="numeric",
                        transformation="list"),
         contains="filter")


## =================================================================
## filterTree
## ----------------------------------------------------------------
## A collection of Filters organized as a DAG. The DAG specifies the
## Filters and the Populations. The Populations are names only because
## they get created only when a FilterTree is applied to a dataset.
## -----------------------------------------------------------------
setClass("filterTree",
         representation(tree="graphNEL",
                        filterSet="list"),
         validity=function(object){
           msg = TRUE
           if(length(edges(tree)) != length(filterSet))
             msg = "\nOne filter must be supplied for each edge."
           return(msg)
         }
         )


## ===========================================================================
## Decision Tree - listed in the standard
## ---------------------------------------------------------------------------
## An object describing a selection applied to a data matrix. Consist of
## a functions that return logical vectors subsetting the data
## ---------------------------------------------------------------------------


## ===========================================================================
## Boolean gate ~ gateSet
## ---------------------------------------------------------------------------
## 
##
## ---------------------------------------------------------------------------



## ===========================================================================
## filterResult
## ---------------------------------------------------------------------------
## A container for the results after applying a filter (e.g. gate) 
## to flow cytometry data with slots 
##  
## 
##   
## 
## ---------------------------------------------------------------------------
setClass("filterResult",
         ## this holds the result of applying a filter to
         ## a single FCS object
         ## the subset is the filter applied to the original data
         ## filterDetails is any output calculated by an algorithmic
         representation(subSet="numeric",
                        filterDetails="list")
         )
