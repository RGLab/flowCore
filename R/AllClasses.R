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
                        parameters="AnnotatedDataFrame",
                        description="list"),
         prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
           parameters=new("AnnotatedDataFrame",
             data=data.frame(name=I(character(0))),
             varMetadata=data.frame(labelDescription="Name in frame",row.names="name")),
           description=list(note="empty")),
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@exprs))
               msg <- "\nslot 'exprs' must be matrix"
             if(!is.list(object@description))
               msg <- "\nslot 'description' must be a list"
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
                   data=data.frame(),
                   varMetadata=data.frame()),
                 colnames=character(0)),
  validity=function(object){
    nc <- length(colnames(object))
	# Why do we do this? It seems like new() should already check this against the representation().
	# replicating it here seems superfluorous.
#	if(!(is(object@phenoData, "AnnotatedDataFrame") && 
#		is(object@colnames, "character") &&
#		is(object@frames, "environment"))) {
#			warning("An element of this object is of the wrong type.")
#			return(FALSE)
#		}

	#Make sure that all of our samples list
	name.check = is.na(match(sampleNames(object),ls(object@frames,all.names=TRUE)))
	if(any(name.check)) {
		name.list = paste(sampleNames(object)[name.check],sep=",")
		warning(paste("These objects are not in the data environment:",name.list))
		return(FALSE)
	}
	
	#Ensure that all frames match our colnames
	if(!all(sapply(sampleNames(object),function(i) {
		x = get(i,env=object@frames)
		if(identical(object@colnames, colnames(x)))  TRUE else { 
			warning(paste(i,"failing colnames check: ",paste(object@colnames,sep=","),"vs",paste(colnames(x),sep=",")))
			FALSE 
		}
	}))) {
		warning("Some items identified in the data environment either have the wrong dimension or type.")
		return(FALSE)
	}
	return(TRUE)
  })





## ===========================================================================
## Virtual filter/subsetting?
## ---------------------------------------------------------------------------
## An object describing a selection applied to a data matrix. Consist of
## a functions that return logical vectors subsetting the data
## ---------------------------------------------------------------------------
setClass("filter", 
         representation("VIRTUAL",filterId="character",
                        parameters="ANY"),
         validity=function(object){
             msg <- TRUE
             if(!is.character(object@filterId) ||
                length(object@filterId)!=1)
               msg <- "\nslot 'filterId' must be character vector of length 1"
             if(!is.character(object@parameters))
             msg <- "\nslot 'parameters' must be  vector"
# Not sure what this test was supposed to do, but it's really annoying
#             test <- matrix(1:length(object@parameters),
#                            ncol=length(object@parameters))
#             colnames(test) <- object@parameters
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
          representation(min="numeric",
                        max="numeric"),
         contains="filter",
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
             if(!is.matrix(object@boundaries) || nrow(object@boundaries)<2)
               msg <- "\nslot 'boundaries' must be character vector longer than 2"
             return(msg)
         })


## ===========================================================================
## Polytope gate
## ---------------------------------------------------------------------------
setClass("polytopeGate",
         contains="filter",
         representation(boundaries="matrix"),
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
setClass("ellipsoidGate",
         representation(focus="matrix",
                        distance="numeric"),
         contains="filter",
         prototype=list(filterId="ALL", focus=matrix(ncol=2, nrow=2),distance=1),
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

## mode gates are rectangle gates that automatically select their dominant region.
setClass("modeGate",representation(bw="ANY",n="numeric"),prototype=list(bw="nrd0",n=512))

## ===========================================================================
## norm2Filter (adapted from prada)
## ---------------------------------------------------------------------------
setClass("norm2Filter",
         ## the slot method holds the method argument to fitNorm2
         ## the slot scale.factor holds the scalefac argument to fitNorm2
         ## transformation holds a list of length giving transformations, if applicable that are
         ## applied to the data before gating
         representation(method="character",
                        scale.factor="numeric",
                        transformation="list",
						n="numeric"),
         contains="filter")

## ===========================================================================
## kmeansFilter 
## ---------------------------------------------------------------------------
setClass("kmeansFilter",
	representation("filter",populations="character"))

## ===========================================================================
## sampleFilter 
## ---------------------------------------------------------------------------
setClass("sampleFilter",representation("filter",size="numeric"))

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

setClass("setOperationFilter",representation("filter",filters="list"))
setClass("unionFilter",representation("setOperationFilter"))
setClass("intersectFilter",representation("setOperationFilter"))
setClass("complementFilter",representation("setOperationFilter"),
	validity=function(object) { 
		if(length(object@filters) != 1) {
			warning("Complement filters can only operate on a single filter")
			return(FALSE)
		}
		TRUE
	})
setClass("subsetFilter",representation("setOperationFilter"),
	validity=function(object) {
		if(length(object@filters) != 2) {
			warning("Subset filters are only defined as binary operators")
			return(FALSE)
		}
		TRUE
	})

## ===========================================================================
## filterResult
## ---------------------------------------------------------------------------
## A container for the results after applying a filter (e.g. gate) 
## to flow cytometry data with slots 
## ---------------------------------------------------------------------------
setClass("filterResult",
	representation("filter",frameId="character",filterDetails="list"),
	prototype=list(frameId=character(0),filterDetails=list()))
setClass("logicalFilterResult",representation("filterResult",subSet="logical"))
setClass("multipleFilterResult",representation("filterResult",subSet="factor"))
setClass("randomFilterResult",representation("filterResult",subSet="numeric"))


## Describing transforms

# Parameterize transforms so that we can describe them.
setClass("transform",representation("function"))

# I want to be able to include transforms within a filter. First we need to know which parameters should
# be input filters
setClass("transformMap",representation(output="character",input="character",f="function"))
#A list of transformMaps
setClass("transformList",representation(transforms="list"))
setClass("transformFilter",representation("filter",transforms="transformList",filter="filter"))


