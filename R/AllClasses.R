## =========================================================================##
## =========================================================================##
##              Class definitions and contructors if available              ##
## =========================================================================##
## =========================================================================##




## ===========================================================================
##  flowFrame
## ---------------------------------------------------------------------------
## A container for flow cytometry measurements with slots exprs, parameters
## and description. exprs contains measurement values, description contains 
## information from file headers of FCS file and parameters contains
## information about the FCS measurement parameters (i.e. channels) available
## ---------------------------------------------------------------------------
setClass("flowFrame",                
         representation(exprs="matrix",
                        parameters="AnnotatedDataFrame",
                        description="list"),
         prototype=list(exprs=matrix(numeric(0), nrow=0, ncol=0),
         parameters=new("AnnotatedDataFrame",
         data=data.frame(name=I(character(0))),
         varMetadata=data.frame(labelDescription="Name in frame",
         row.names="name")),
         description=list(note="empty")))



## ===========================================================================
##  flowSet
## ---------------------------------------------------------------------------
## A collection of several cytoFrames making up one experiment. Slots 
## frames, phenoData, colnames. Frames contains the cytoFrame objects,
## phenoData the experiment meta data and colnames the channel names. 
## ---------------------------------------------------------------------------
setClass("flowSet",                   
         representation(frames="environment",
                        phenoData="AnnotatedDataFrame",
                        colnames="character"),
         prototype=list(frames=new.env(hash=TRUE, parent=emptyenv()),
         phenoData=new("AnnotatedDataFrame",
         data=data.frame(),
         varMetadata=data.frame()),
         colnames=character(0)),
         validity=function(object){
             nc <- length(colnames(object))
             ## Make sure that all of our samples list
             name.check <- is.na(match(sampleNames(object),ls(object@frames,
                                                              all.names=TRUE)))
             if(any(name.check)) {
                 name.list <- paste(sampleNames(object)[name.check],sep=",")
                 return(paste("These objects are not in the data environment:",name.list))
             }
             ##Ensure that all frames match our colnames
             if(!all(sapply(sampleNames(object),function(i) {
                 x = get(i,env=object@frames)
                 if(identical(object@colnames, colnames(x))){
                     TRUE
                 }else{ 
                     return(paste(i, "failing colnames check: ",
                                  paste(object@colnames, sep=","),
                                  "vs", paste(colnames(x), sep=",")))
                 }
             }))){
                 return("Some items identified in the data environment either ",
                        "have the wrong dimension or type.")
             }
             return(TRUE)
         })

## constructor
flowSet = function(...,phenoData) {
    x = list(...)
    if(length(x) == 1 && is.list(x[[1]])) x = x[[1]]
    f = as(x,"flowSet")
    if(!missing(phenoData)) phenoData(f) = phenoData
    f
}



## ===========================================================================
## Virtual filter and derived concreteFilter and parameterFilter
## ---------------------------------------------------------------------------
## A class describing a selection applied to a data matrix. Consist of
## a filterId and the names of the parameters to operate on (for parameter
## filters only). Specific filters all inherit from either of these two
## classes
## ---------------------------------------------------------------------------
setClass("filter", 
         representation("VIRTUAL", filterId="character"),
         prototype=prototype(filterId=""),
         validity=function(object) {
             if(length(object@filterId) != 1)
                 "'filterId' must have a length of one."
             else
                 TRUE
         })
setClass("concreteFilter", "filter")
setClass("parameterFilter", representation(parameters="character"),
         contains="concreteFilter", prototype=prototype(parameters=""))



## ===========================================================================
## Rectangular gate
## ---------------------------------------------------------------------------
## A class describing a 2D rectangular region in the parameter space. Slots
## min and max hold the boundaries in the two dimensions
## ---------------------------------------------------------------------------
setClass("rectangleGate",
         representation(min="numeric",
                        max="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="Rectangle Gate",
         min=0,max=Inf)
         )

## constructor
rectangleGate <- function(filterId="rectangleGate", .gate,...) {
    if(missing(.gate) || !is.matrix(.gate))
      	.gate <- sapply(if(missing(.gate)) list(...) else .gate,function(x) {
            x = sort(x);c("min"=x[1],"max"=x[2])
		})
    new("rectangleGate", filterId=filterId, parameters=colnames(.gate),
        min=.gate[1,], max=.gate[2,])
}



## ===========================================================================
## Polygon gate
## ---------------------------------------------------------------------------
## A class describing a 2D polygonal region in the parameter space. Slot
## boundary holds the vertices of the polygon in a 2 colum matrix
## ---------------------------------------------------------------------------
setClass("polygonGate",
         representation(boundaries="matrix"),
         contains="parameterFilter",
         prototype=list(filterId="ALL", boundaries=matrix(ncol=2, nrow=3)),
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@boundaries) || nrow(object@boundaries)<2)
                 msg <- "\nslot 'boundaries' must be a numeric matrix of at least 2 rows"
             return(msg)
         })

## constructor
polygonGate <- function(filterId="polygonGate", boundaries,...) {
    if(missing(boundaries) || !is.matrix(boundaries)) 
        boundaries = {as.matrix(if(missing(boundaries)) do.call("cbind",list(...))
        else boundaries)}
    new("polygonGate",filterId=filterId, parameters=colnames(boundaries),
        boundaries=boundaries)
}



## ===========================================================================
## Polytope gate
## ---------------------------------------------------------------------------
## A class describing a 2D polytope region in the parameter space. Slot
## boundary holds the vertices of the polygon in a 2 colum matrix
## ---------------------------------------------------------------------------
setClass("polytopeGate",
         representation(boundaries="matrix"),
         contains="parameterFilter",
         prototype=list(filterId="ALL", boundaries=matrix(ncol=2, nrow=2)), 
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@boundaries) ||
                nrow(object@boundaries)<2)
               msg <- "\nslot 'boundaries' must be a numeric matrix of at least 2 rows"
             return(msg)
         })

## constructor
polytopeGate <- function(filterId="polytopeGate", .gate, ...) {
    if(missing(.gate) || !is.matrix(.gate))
        .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
    new("polytopeGate", filterId=filterId, parameters=colnames(.gate),
        boundaries=.gate)
}



## ===========================================================================
## Ellipsoid gate
## ---------------------------------------------------------------------------
## A class describing a 2D polytope region in the parameter space. Slots
## focus and distance hold the xy position of the centroid and the distance
## respectively.
## ---------------------------------------------------------------------------
setClass("ellipsoidGate",
         representation(focus="matrix",
                        distance="numeric"),
         contains="parameterFilter",
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

ellipsoidGate <- function(filterId="ellipsoidGate", .gate, distance,...) {
    if(missing(.gate) || !is.matrix(.gate))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
      
    new("ellipsoidGate", filterId=filterId, parameters=colnames(.gate),
        focus=.gate, distance=distance)
}



## ===========================================================================
## norm2Filter
## ---------------------------------------------------------------------------
## A class to describe the fit of a bivariate normal distribution.
## Slot method holds the method argument to fitNorm2,
## slot scale.factor holds the scalefac argument to fitNorm2,
## transformation holds a list of length giving transformations, if applicable
## that are applied to the data before gating. n is the number of points used
## in the subsampling step.
## ---------------------------------------------------------------------------
setClass("norm2Filter",
         representation(method="character",
                        scale.factor="numeric",
                        transformation="list",
                        n="numeric"),
         contains="parameterFilter")

## constructor
norm2Filter <- function(x, y, method="covMcd", scale.factor=1,
                        filterId="norm2Gate", n=50000,...)
{
    if(missing(y)) {
        if(length(x)==1)
            stop("You must specify two parameters for a norm2 gate.")
        if(length(x)>2)
            warning("Only the first two parameters will be used.")
        y=x[2]
        x=x[1]
    } else {
        if(length(x)>1 || length(y)>1)
            warning("Only the first two parameters from 'x' and ",
                    "'y' will be used.")
        x = x[1]
        y = y[1]
    }
    new("norm2Filter",parameters=c(x,y),method=method,scale.factor=scale.factor,
        filterId=filterId,n=n,...)
}


## ===========================================================================
## kmeansFilter 
## ---------------------------------------------------------------------------
setClass("kmeansFilter",
         representation(populations="character"),
         contains="parameterFilter")

## contructor
## not sure why but the parameters in list format can not be read.
kmeansFilter = function(filterId="kmeans",...)
{
    l = length(list(...))
    if(l>1)
	stop("k-means filters only operate on a single parameter.")
    x = ..1
    if(is.list(x)) {
        new("kmeansFilter",parameters=names(x)[1],populations=x[[1]],
            filterId=filterId)
    } else {
        new("kmeansFilter",parameters=names(list(...))[1],populations=x,
            filterId=filterId)
    }
}



## ===========================================================================
## curv2Filter
## this filter can hold parameters to find siginficant high density regions
## in two dimensions based on Matt Wand's feature software
## ---------------------------------------------------------------------------
setClass("curv2Filter",
         representation(bwFac="numeric",
                        gridsize="numeric"),
         contains="parameterFilter")

##constructor
curv2Filter <-
    function(x, y, filterId="curv2Filter", bwFac=1.2,
             gridsize=rep(151,2), ...)
{
    if(!is.numeric(bwFac) || length(bwFac)!=1)
        stop("'bwFac must be numeric skalar")
    if(!is.numeric(gridsize) || length(gridsize)!=2)
        stop("'gridsize must be numeric skalar")
    if(missing(y)) {
        if(length(x)==1)
            stop("You must specify two parameters for a curv2 filter.")
        if(length(x)>2)
            warning("Only the first two parameters will be used.")
        y=x[2]
        x=x[1]
    } else {
        if(length(x)>1 || length(y)>1)
            warning("Only the first two parameters from 'x' and ",
                    "'y' will be used.")
        x = x[1]
        y = y[1]
    }
    new("curv2Filter",parameters=c(x,y), bwFac=bwFac, gridsize=gridsize,
        filterId=filterId, ...)
}

   

## ===========================================================================
## sampleFilter 
## ---------------------------------------------------------------------------
setClass("sampleFilter",
         representation(size="numeric"),
         contains="concreteFilter")

sampleFilter = function(filterId="sample",size) new("sampleFilter",filterId=filterId,size=size)


## ===========================================================================
## expressionFilter 
## Let's us encapsulate an expression as a gate
## ---------------------------------------------------------------------------
setClass("expressionFilter",
         representation(expr="call",args="list"),
         contains="concreteFilter")
expressionFilter = function(expr,...,filterId) {
	expr = substitute(expr)
	if(missing(filterId)) filterId = deparse(expr)
	new("expressionFilter",filterId=filterId,expr=expr,args=list(...))
}


## ===========================================================================
## filterSet 
## ---------------------------------------------------------------------------
setClass("filterSet",representation(env="environment"),prototype=prototype(env=new.env(hash=TRUE,parent=emptyenv())))
filterSet = function(...) {
	filters = list(...)
	#Allow the list(x,y,z) format as well.
	if(length(filters)==1 && is.list(filters[[1]])) filters = filters[[1]]
	if(length(filters) == 0)
		new("filterSet",env=new.env(parent=emptyenv()))
	else
		as(filters,"filterSet")
}
#References a filter (contained within a filterSet)
setClass("filterReference",representation(name="character",env="environment"),contains="filter")

## ===========================================================================
## setOperationFilter
## ---------------------------------------------------------------------------
setClass("setOperationFilter",
         representation(filters="list"),
         contains="concreteFilter")


## ===========================================================================
## unionFilter 
## ---------------------------------------------------------------------------
setClass("unionFilter",representation("setOperationFilter"))


## ===========================================================================
## intersectFilter 
## ---------------------------------------------------------------------------
setClass("intersectFilter",representation("setOperationFilter"))


## ===========================================================================
## complementFilter 
## ---------------------------------------------------------------------------
setClass("complementFilter",representation("setOperationFilter"),
	validity=function(object) { 
		if(length(object@filters) != 1) {
			warning("Complement filters can only operate on a single filter")
			return(FALSE)
		}
		TRUE
	})


## ===========================================================================
## subsetFilter 
## ---------------------------------------------------------------------------
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
         representation(frameId="character", filterDetails="list"),
         contains="concreteFilter",
         prototype=list(frameId=character(0), filterDetails=list()))


## ===========================================================================
## logicalFilterResult
## ---------------------------------------------------------------------------
setClass("logicalFilterResult",
         representation(subSet="logical"),
         contains="filterResult")


## ===========================================================================
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("multipleFilterResult",
         representation(subSet="factor"),
         contains="filterResult")


## ===========================================================================
## manyFilterResult
## ---------------------------------------------------------------------------
## A special case of multipleFilterResult that arises when there are
## overlapping sets
## ---------------------------------------------------------------------------
setClass("manyFilterResult",
         representation(subSet="matrix",dependency="ANY"),
         contains="filterResult")
manyFilterResult = function(filters,frameId,dependency=NULL) {
	q = new("manyFilterResult",
		filterDetails = sapply(filters,slot,"filterDetails"),
		subSet=do.call("cbind",lapply(filters,as,"logical")),
		dependency=dependency)
	colnames(q@subSet) = sapply(filters,slot,"filterId")
	q
}
## ===========================================================================
## randomFilterResult
## ---------------------------------------------------------------------------
setClass("randomFilterResult",
         representation(subSet="numeric"),
         contains="filterResult")

#
# filterSummary now becomes a legitimate class.
setClass("filterSummary",representation(name="character",true="numeric",count="numeric",p="numeric"))



## ===========================================================================
## transform
## ---------------------------------------------------------------------------
## Parameterize transforms so that we can describe them.
## ---------------------------------------------------------------------------
setClass("transform", representation(transformationId="character", "function"))

## Linear transformation function
linearTransform <- function(transformationId,a=1,b=0){
    if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
    t= new("transform", .Data=function(x){    
        x <- a*x+b
    })
    t@transformationId = transformationId
    t
  }

## Quadratic transformation function
quadraticTransform <- function(transformationId,a=1,b=1,c=0){
  if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
  if(!is.double(c))
       stop("c must be numeric")
    t = new("transform",.Data=function(x){
        x <- a*x^2 + b*x + c
      })
  t@transformationId = transformationId
  t
}

## Natural logarithm transformation function
lnTransform <- function(transformationId,r=1,d=1){
    if(!is.double(r) || r <= 0)
       stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
       stop("d must be numeric")
    t= new("transform",.Data=function(x){
      x<-log(x)*(r/d)
    })
    t@transformationId = transformationId
    t
}

## Logarithm transformation function
logTransform <- function(transformationId,logbase=10,r=1,d=1){

  if(!is.double(r) || r <= 0)
    stop("r must be numeric and positive")
  if(!is.double(d) || d <=0)
    stop("d must be numeric")
  if(!is.double(r) || r <=0)
    stop("r must be numeric and positive")
  if(!is.double(logbase) || logbase <= 1)
    stop("logabse must be a pnumeric greater than 1")
  t = new("transform",.Data=function(x){
    x <- log(x,logbase)*(r/d)
  })
  t@transformationId = transformationId
  t
}


## General biexponential transformation function
biexponentialTransform<- function(transformationId, a=.5, b=1,c=.5,d=1,f=0,w=0,
           tol=.Machine$double.eps^0.25,maxit=as.integer(5000)){
    
    t = new("transform",.Data=function(x){
        x <- .Call(biexponential_transform, x, a, b, c, d, f, w, tol, maxit)
    })
    t@transformationId = transformationId
    t
}

## Logicle transformation function
logicleTransform <- function(transformationId, w=0,r=262144,d=5,...) {
  if(w>d) stop("Negative range decades must be smaller than total number of decades")
    w = w*log(10)
    d = d*log(10)
    if(w==0) p = 1 else p = uniroot(function(p) -w+2*p*log(p)/(p+1), c(.Machine$double.eps,2*(w+d)))$root
    t= new("transform",.Data=biexponentialTransform(transformationId, a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...))
    t@transformationId = transformationId
    t
  }

## Truncation Transformation
truncateTransform <- function(transformationId, a=1){
  t= new("transform",.Data=function(x){
    x[x<=a] <- a
    x
  })
  t@transformationId = transformationId
  t
}

## Scale Transformation
scaleTransform <- function(transformationId,a=1,b=10^4){
    t = new("transform",.Data=function(x){
     	x=(x-a)/(b-a)
    })
    t@transformationId = transformationId
    t
}

## Split-scale Transformation
splitScaleTransform <- function(transformationId, maxValue=1023,
                                transitionChannel=64, r=192){
    
    maxChannel = r + transitionChannel
    b = transitionChannel/2
    d= 2*log10(exp(1))*r/transitionChannel
    logt = -2*log10(exp(1))*r/transitionChannel + log10(maxValue)
    t = 10^logt
    a = transitionChannel/(2*t)
    logCT = (a*t+b)*d/r
    c = 10^logCT/t
    
    tr = new("transform",.Data= function(x){
        idx = which(x <= t)
        idx2 = which(x > t)

        if(length(idx2)>0){
            x[idx2]= log10(c*x[idx2])*r/d
        }
        if(length(idx)>0){
            x[idx] = a*x[idx]+b
        }
        x        
    })
    tr@transformationId = transformationId
    tr
}

## Hyperbolic Arcsin Transformation
arcsinhTransform <- function(transformationId,a=1,b=1,c=0) {
  t = new("transform",.Data=function(x) asinh(a+b*x)+c)
  t@transformationId = transformationId
  t
}

## ===========================================================================
## parameterTransform
## ---------------------------------------------------------------------------
## A class that only applied a parameter transform to a subset of the
## parameters during an %on% operation.
## ---------------------------------------------------------------------------
setClass("parameterTransform",representation(parameters="character"),contains="transform")
parameterTransform = function(FUN,params)
	new("parameterTransform",.Data=as.function(FUN),parameters=as.character(params))

## ===========================================================================
## transformMap
## ---------------------------------------------------------------------------
## I want to be able to include transforms within a filter. First we need to
## know which parameters should be input filters
## ---------------------------------------------------------------------------
setClass("transformMap",
         representation(output="character", input="character", f="function"))


## ===========================================================================
## transformList
## ---------------------------------------------------------------------------
## A list of transformMaps
## ---------------------------------------------------------------------------
setClass("transformList", representation(transforms="list"))


## ===========================================================================
## transformFilter
## ---------------------------------------------------------------------------
setClass("transformFilter",
         representation(transforms="transformList", filter="filter"),
         contains="concreteFilter")

