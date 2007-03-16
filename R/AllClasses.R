## Class definitions and contructors if available


## ===========================================================================
##  flowFrame
## ---------------------------------------------------------------------------
## A container for flow cytometry measurements with slots exprs, parameters
## and description. exprs contains measurement values, description contains 
## information from file headers of FCS file and parameters contains
## information about the different FCS measurement parameters (i.e. channels)
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
           description=list(note="empty")),
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@exprs))
               msg <- "\nslot 'exprs' must be matrix"
             if(!is.list(object@description))
               msg <- "\nslot 'description' must be a list"
         })


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
         prototype=list(frames=new.env(),
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
             warning(paste("These objects are not in the data environment:",
                           name.list))
             return(FALSE)
           }
           ##Ensure that all frames match our colnames
           if(!all(sapply(sampleNames(object),function(i) {
             x = get(i,env=object@frames)
             if(identical(object@colnames, colnames(x))){
               TRUE
             }else{ 
               warning(paste(i, "failing colnames check: ",
                             paste(object@colnames, sep=","),
                             "vs", paste(colnames(x), sep=",")))
               FALSE 
             }
           }))){
             warning("Some items identified in the data environment either ",
                     "have the wrong dimension or type.")
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
             if(!is.character(object@parameters))#why is the sig for parameters "ANY" but we test for character?
             msg <- "\nslot 'parameters' must be a character vector"
            return(msg)
         })


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

rectangleGate <- function(filterId="rectangleGate", .gate,...) {
    if(missing(.gate) || !is.matrix(.gate))
      	.gate <- sapply(if(missing(.gate)) list(...) else .gate,function(x)
                        c("min"=x[1],"max"=x[2]))
	new("rectangleGate",filterId=filterId,parameters=colnames(.gate),
            min=.gate[1,],max=.gate[2,])
}


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
               msg <- "\nslot 'boundaries' must be a numeric matrix of at least 2 rows"
             return(msg)
         })

polygonGate <- function(filterId="polygonGate", boundaries,...) {
    if(missing(boundaries)){
      boundaries <- matrix(ncol=2, nrow=3)
      colnames(boundaries) <- rep(NA,2)
    }
    if(!is.matrix(boundaries))
        boundaries <- sapply(if(missing(boundaries)) list(...) else boundaries,
                           function(x) x)
    new("polygonGate",filterId=filterId, parameters=colnames(boundaries),
        boundaries=boundaries)
}


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
               msg <- "\nslot 'boundaries' must be a numeric matrix of at least 2 rows"
             return(msg)
         })

polytopeGate <- function(filterId="polytopeGate", .gate, ...) {
    if(missing(.gate) || !is.matrix(.gate))
      ##nrowGate <- max(unlist(lapply(list(...),length)))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
         
    new("polytopeGate", filterId=filterId, parameters=colnames(.gate),
        boundaries=.gate)
}


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

ellipsoidGate <- function(filterId="ellipsoidGate", .gate, distance,...) {
    if(missing(.gate) || !is.matrix(.gate))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
      
    new("ellipsoidGate", filterId=filterId, parameters=colnames(.gate),
        focus=.gate, distance=distance)
}

## ===========================================================================
## norm2Filter
## ---------------------------------------------------------------------------
## the slot method holds the method argument to fitNorm2
## the slot scale.factor holds the scalefac argument to fitNorm2
## transformation holds a list of length giving transformations, if applicable
## that are applied to the data before gating
## ---------------------------------------------------------------------------
setClass("norm2Filter",
         representation(method="character",
                        scale.factor="numeric",
                        transformation="list",
                        n="numeric"),
         contains="filter")

norm2Filter <- function(x,y,method="covMcd",scale.factor=1,filterId="norm2Gate",
                        n=50000,...) {
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
            filterId=filterId,n=50000,...)
}


## ===========================================================================
## kmeansFilter 
## ---------------------------------------------------------------------------
setClass("kmeansFilter",
         representation(populations="character"),
         contains="filter")

kmeansFilter = function(filterId="kmeans",...) {
	l = length(list(...))
	if(l>1)
		stop("k-means filters only operate on a single parameter.")
	x = ..1
	if(is.list(x)) {
		new("kmeansFilter",parameters=names(x)[1],populations=x[[1]],filterId=filterId)
	} else
		new("kmeansFilter",parameters=names(list(...))[1],populations=x,filterId=filterId)
}


## ===========================================================================
## sampleFilter 
## ---------------------------------------------------------------------------
setClass("sampleFilter",
         representation(size="numeric"),
         contains="filter")

sampleFilter = function(filterId="sample",size) {
	new("sampleFilter",parameters=character(0),filterId=filterId,size=size)
}


## ===========================================================================
## multiFilter 
## ---------------------------------------------------------------------------
setClass("multiFilter",
         representation(populations="character",filters="list"),
         contains="filter")



## =================================================================
## filterCollection
## ----------------------------------------------------------------
## A collection of Filters organized as a DAG. The DAG specifies the
## Filters and the Populations. The Populations are names only because
## they get created only when a FilterTree is applied to a dataset.
## -----------------------------------------------------------------
## Environment can't be extended
## setClass("filterCollection",representation("environment"))

setClass("filterCollection",representation("list"))
## ===========================================================================
## promisedFilter 
## ---------------------------------------------------------------------------
## A reference to a filter in an environment.
## ---------------------------------------------------------------------------
setClass("promisedFilter",representation("filter",name="character"))



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
## setOperationFilter
## ---------------------------------------------------------------------------
setClass("setOperationFilter",
         representation(filters="list"),
         contains="filter")


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
         contains="filter",
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
         representation(subSet="logical"),
         contains="filterResult")

## ===========================================================================
## randomFilterResult
## ---------------------------------------------------------------------------
setClass("randomFilterResult",
         representation(subSet="numeric"),
         contains="filterResult")


## ===========================================================================
## transform
## ---------------------------------------------------------------------------
## Parameterize transforms so that we can describe them.
## ---------------------------------------------------------------------------
setClass("transform", representation("function"))

## Linear transformation function
linearTransform <- function(transformationId,a=1,b=0){
    if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
    new("transform",.Data=function(x){    
        x <- a*x+b
    })
}

## Quadratic transformation function
quadraticTransform <- function(transformationId,a=1,b=1,c=0){
  if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
  if(!is.double(c))
       stop("c must be numeric")
    new("transform",.Data=function(x){
        x <- a*x^2 + b*x + c
    })
}

## Natural logarithm transformation function
lnTransform <- function(transformationId,r=1,d=1){
    if(!is.double(r) || r <= 0)
       stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
       stop("d must be numeric")
    new("transform",.Data=function(x){
     x<-log(x)*(r/d)
 	})
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
    new("transform",.Data=function(x){
        x <- log(x,logbase)*(r/d)
    })
}


## General biexponential transformation function
biexponentialTransform<- function(transformationId, a=.5, b=1,c=.5,d=1,f=0,w=0,
           tol=.Machine$double.eps^0.25,maxit=as.integer(5000)){
    new("transform",.Data=function(x){
        x <- .Call(biexponential_transform,x,a,b,c,d,f,w,tol,maxit)
    })
}

## Logicle transformation function
logicleTransform <- function(transformationId, w=0,r=262144,d=5,...) {
    if(w>d) stop("Negative range decades must be smaller than total number of decades")
    w = w*log(10)
    d = d*log(10)
    if(w==0) p = 1 else uniroot(function(p) -w+2*p*log(p)/(p+1), c(.Machine$double.eps,2*(w+d)))$root
    new("transform",.Data=function(x)
        biexponentialTransform(transformationId, a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...))
}
## Truncation Transformation
truncateTransform <- function(transformationId,a=1){
    new("transform",.Data=function(x){
        x[x<a] <- a
        x
    })
}

## Scale Transformation
scaleTransform <- function(transformationId,a=1,b=10^4){
    new("transform",.Data=function(x){
     	x=(x-a)/(b-a)
    })
}

## Hyperbolic Arcsin Transformation
arcsinhTransform <- function(transformationId,a=1,b=1,c=0) {
	new("transform",.Data=function(x) asinh(a+b*x)+log(c))
}



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
         contains="filter")


