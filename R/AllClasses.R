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
         parameters=new("AnnotatedDataFrame"),
         description=list(note="empty")))


## helper function to create empty AnnotatedDataFrame for the parameters slot
parDefault <- function(exp){
    vm <- data.frame(labelDescription=c(name="Name of Parameter",
                     desc="Description of Parameter",
                     range="Range of Parameter",
                     minRange="Minimum Parameter Value after Transformation",
                     maxRange="Maximum Parameter Value after Transformation"))
    pd <- data.frame(name=colnames(exp), desc=colnames(exp),
                     range=apply(exp, 2, max, na.rm=TRUE),
                     minRange=apply(exp, 2, min, na.rm=TRUE),
                     maxRange=apply(exp, 2, max, na.rm=TRUE))
    new("AnnotatedDataFrame", pd, vm)
}

## check parameter AnnotatedDataFrame for validity
isValidParameters <- function(parms, exprs)
{
     if(!is(parms, "AnnotatedDataFrame"))
           stop("Argument 'parameters' must be object of class ",
                'AnnotatedDataFrame', call.=FALSE)
     if(!all(c("name", "desc", "range", "minRange", "maxRange")
             %in% varLabels(parms)))
         stop("The following columns are mandatory:\n  'name', 'desc',",
              "'range', 'minRange', 'maxRange'", call.=FALSE)
     if(!missing(exprs))
         if(!all(colnames(exprs) %in% parms$name))
             stop("parameter description doesn't match colnames of the ",
                  "data matrix", call.=FALSE)
     return(TRUE)
}

## constructor
flowFrame <- function(exprs, parameters, description=list()){
    if(!is.matrix(exprs) || !is.numeric(exprs) || is.null(colnames(exprs)))
        stop("Argument 'exprs' must be numeric matrix with colnames ",
             "attribute set", call.=FALSE)
    if(missing(parameters))
        parameters <- parDefault(exprs)
    else
        isValidParameters(parameters, exprs)
    if(!is.list(description))
        stop("Argument 'description' must be a list", call.=FALSE)
    new("flowFrame", exprs=exprs, parameters=parameters,
        description=description)
}
  

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
             name.check <- is.na(match(sampleNames(object), ls(object@frames,
                                                              all.names=TRUE)))
             if(any(name.check)) {
                 name.list <- paste(sampleNames(object)[name.check], sep=",")
                 return(paste("These objects are not in the data environment:",
                              name.list))
             }
             ##Ensure that all frames match our colnames
             if(!all(sapply(sampleNames(object), function(i) {
                 x <- get(i, env=object@frames)
                 if(all(object@colnames %in% colnames(x))){
                     TRUE
                 }else{ 
                     return(paste(i, "failing colnames check: ",
                                  paste(object@colnames, sep=","),
                                  "vs", paste(colnames(x), sep=",")))
                 }
             }))){
                 return(paste("Some items identified in the data environment",
                              "either have the wrong dimension or type."))
             }
             return(TRUE)
         })

## constructor
flowSet <- function(..., phenoData) {
    x <- list(...)
    if(length(x) == 1 && is.list(x[[1]]))
        x <- x[[1]]
    if(!all(sapply(x, is, "flowFrame")))
        stop("All additional arguments must be flowFrames")
    f <- as(x, "flowSet")
    if(!missing(phenoData))
        phenoData(f) = phenoData
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

setClass("concreteFilter",
         contains="filter")

setClass("parameterFilter",
         representation(parameters="character"),
         contains="concreteFilter",
         prototype=prototype(parameters=""))



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
         min=0, max=Inf)
         )

## constructor
rectangleGate <- function(filterId="rectangleGate", .gate, ...) {
    if(missing(.gate) || !is.matrix(.gate))
      	.gate <- sapply(if(missing(.gate)) list(...) else .gate,
                        function(x){
                            x <- sort(x)
                            c("min"=x[1], "max"=x[2])
                        })
    if(!is.numeric(.gate) || is.null(colnames(.gate)))
        stop("Expecting named arguments or a single named matrix\n",
             "as input for gate boundaries.", call.=FALSE)
    new("rectangleGate", filterId=filterId, parameters=colnames(.gate),
        min=.gate[1,], max=.gate[2,])
}




## ===========================================================================
## Quadrant gate
## ---------------------------------------------------------------------------
## A class describing a gate which separated a 2D paramerter space into
## four quadrants. Slot boundary holds a vector of length two indicating
## the quadrant boundaries in each of the two dimensions
## ---------------------------------------------------------------------------
setClass("quadGate",
         representation(boundary="numeric"),        
         contains="parameterFilter",
         prototype=list(filterId="Quadrant Gate", c(Inf, Inf)))

## constructor
quadGate <- function(filterId="quadGate", .gate, ...) {
    if(missing(.gate))
        .gate <- list(...)
    n <- names(.gate) 
    .gate <- as.numeric(.gate)
    names(.gate) <- n
    if(!is.numeric(.gate) || length(names(.gate)) !=2 |
       any(names(.gate)==""))
        stop("Expecting two named arguments or a single named vector\n",
             "of length 2 as input for gate boundaries.", call.=FALSE)
    new("quadGate", filterId=filterId, parameters=names(.gate),
        boundary=.gate)
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
             if(!is.matrix(object@boundaries) || nrow(object@boundaries)<3 ||
                ncol(object@boundaries)!=2)
                 msg <- paste("\nslot 'boundaries' must be a numeric matrix",
                              "of at least 3 rows and exactly 2 columns")
             return(msg)
         })

## constructor
polygonGate <- function(filterId="polygonGate", boundaries, ...) {
    if(missing(boundaries) || !is.matrix(boundaries)) 
        boundaries = {as.matrix(if(missing(boundaries))
                                do.call("cbind", list(...))
        else boundaries)}
    new("polygonGate", filterId=filterId, parameters=colnames(boundaries),
        boundaries=boundaries)
}



## ===========================================================================
## Polytope gate
## ---------------------------------------------------------------------------
## A class describing a nD polytope region in the parameter space. Slot
## boundary holds the vertices of the polygon in a n colum matrix
## ---------------------------------------------------------------------------
## FIXME: Need a representation of halfspaces and corresponding linear
## equations here.

setClass("polytopeGate",
         representation(boundaries="matrix"),
         contains="parameterFilter",
         prototype=list(filterId="ALL", boundaries=matrix()), 
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@boundaries) || nrow(object@boundaries)<3)
                 msg <- paste("\nslot 'boundaries' must be a numeric matrix",
                              "of at least 2 rows")
             return(msg)
         })

## constructor
polytopeGate <- function(filterId="polytopeGate", .gate, ...) {
    if(missing(.gate) || !is.matrix(.gate))
        .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
    if(is.null(colnames(.gate)))
        stop("either provide colnames for argument '.gate' or use named ",
             "arguments for boundary vertices", call.=FALSE)
    new("polytopeGate", filterId=filterId, parameters=colnames(.gate),
        boundaries=.gate)
}



## ===========================================================================
## Ellipsoid gate
## ---------------------------------------------------------------------------
## A class describing an ellipsoid region in the parameter space. Slots
## mean and cov contain the mean values and the covariance matrix describing
## the ellipse
## ---------------------------------------------------------------------------
setClass("ellipsoidGate",
         representation(mean="numeric",
                        cov="matrix",
			distance="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="ALL", mean=numeric(), cov=matrix(), distance=1),
         validity=function(object){
             msg <- TRUE
             if(!is.matrix(object@cov) ||
                nrow(object@cov) != ncol(object@cov) ||
                nrow(object@cov) < 2) 
                 msg <- "\nslot 'cov' must be a symmetric matrix of at least 2 rows"
             if(!is.numeric(object@mean) ||
                length(object@mean) != nrow(object@cov))
                 msg <- paste("\nslot 'mean' must be numeric vector of",
                              "same length as dimensions in 'cov'")
	     if(!is.numeric(object@distance) ||	length(object@distance)!=1)
	         msg <- "'distance' must be numeric of length 1"      
             return(msg)
         })

## constructor
ellipsoidGate <- function(filterId="ellipsoidGate", .gate, mean, distance=1,
...) {
    if(missing(.gate) || !is.matrix(.gate))
      .gate <- sapply(if(missing(.gate)) list(...) else .gate, function(x) x)
    rn <- rownames(.gate)
    cn <- colnames(.gate)
    if(any(sapply(dimnames(cov), is.null)) ||
       !all(rn == cn) || !all(names(mean) == cn))
        stop("'.gate' must be covariance matrix with same dimnames as ",
             "parameter names in 'mean'", call.=FALSE)
    names(mean) <- cn
    new("ellipsoidGate", filterId=filterId, parameters=colnames(.gate),
        cov=.gate, mean=mean, distance=distance)
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
                        filterId="norm2Gate", n=50000)
{
    x <- unlist(x)
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
    new("norm2Filter", parameters=c(x, y), method=method,
        scale.factor=scale.factor, filterId=filterId, n=n)
}


## ===========================================================================
## kmeansFilter
## ---------------------------------------------------------------------------
## Apply kmeans clustering on a single parameter. The number k of clusters
## is given by the length of the 'populations' slot. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("kmeansFilter",
         representation(populations="character"),
         contains="parameterFilter")

## contructor
kmeansFilter <- function(..., filterId="kmeans")
{
    ll <- list(...)
    l <- length(ll)
    if(l>1)
	warning("k-means filters only operate on a single parameter.\n",
                "Using '", names(ll)[1], "' now.", call.=FALSE)
    x <- ..1
    if(is.character(x))
        new("kmeansFilter", parameters=names(ll)[1], populations=x,
            filterId=filterId)
    else if(is.list(x) && ! length(names(x))){
        new("kmeansFilter", parameters=names(ll)[1], populations=unlist(x),
            filterId=filterId)
    }else{
        if(length(x)>1){
            warning("k-means filters only operate on a single parameter.\n",
                    "Using '", names(x)[1], "' now.", call.=FALSE)
            x <- x[1]
        }
        new("kmeansFilter", parameters=names(x[1]), populations=unlist(x),
            filterId=filterId)
    }
}



## ===========================================================================
## curv2Filter
## ---------------------------------------------------------------------------
## this filter can hold parameters to find siginficant high density regions
## in two dimensions based on Matt Wand's feature software. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("curv2Filter",
         representation(bwFac="numeric",
                        gridsize="numeric"),
         contains="parameterFilter")

##constructor
curv2Filter <-
    function(x, y, filterId="curv2Filter", bwFac=1.2,
             gridsize=rep(151, 2))
{
    if(!is.numeric(bwFac) || length(bwFac)!=1)
        stop("'bwFac must be numeric skalar.", call.=FALSE)
    if(!is.numeric(gridsize) || length(gridsize)!=2)
        stop("'gridsize must be numeric skalar", call.=FALSE)
    if(missing(y)) {
        if(length(x)==1)
            stop("You must specify two parameters for a curv2Filter.",
                 call.=FALSE)
        if(length(x)>2)
            warning("Only using parameters '", x[1], "' and '", x[2],
                    "'.", call.=FALSE)
        y=x[2]
        x=x[1]
    } else {
        if(length(x)>1 || length(y)>1)
            warning("Only using parameters '", x[1], "' and '", y[1],
                    "'.", call.=FALSE)
        x = x[1]
        y = y[1]
    }
    new("curv2Filter", parameters=as.character(c(x, y)), bwFac=bwFac,
        gridsize=gridsize, filterId=as.character(filterId))
}



## ===========================================================================
## curv1Filter
## ---------------------------------------------------------------------------
## this filter can hold parameters to find siginficant high density regions
## in one dimension based on Matt Wand's feature software. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("curv1Filter",
         representation(bwFac="numeric",
                        gridsize="numeric"),
         contains="parameterFilter")

##constructor
curv1Filter <-
    function(x, filterId="curv1Filter", bwFac=1.2,
             gridsize=rep(151, 2))
{
    if(!is.numeric(bwFac) || length(bwFac)!=1)
        stop("'bwFac must be numeric skalar")
    if(!is.numeric(gridsize) || length(gridsize)!=2)
        stop("'gridsize must be numeric skalar")
    x <- unlist(x)
    if(length(x)!=1)
        warning("You must specify a single parameters for a curv1Filter.\n",
                "Only using parameter '", x[1], "' now.", call.=FALSE)
    new("curv1Filter", parameters=as.character(x[1]), bwFac=bwFac,
        gridsize=gridsize, filterId=as.character(filterId))
}

   

## ===========================================================================
## sampleFilter
## ---------------------------------------------------------------------------
## sample 'size' rows from a flowFrame
## ---------------------------------------------------------------------------
setClass("sampleFilter",
         representation(size="numeric"),
         contains="concreteFilter")

##constructor
sampleFilter <- function(filterId="sample", size) new("sampleFilter",
                        filterId=filterId, size=size)



## ===========================================================================
## expressionFilter
## ---------------------------------------------------------------------------
## Let's us encapsulate an expression as a gate
## ---------------------------------------------------------------------------
setClass("expressionFilter",
         representation(expr="expression", args="list", deparse="character"),
         contains="concreteFilter")

## constructor
expressionFilter <- function(expr, ..., filterId)
{
    subs <- substitute(expr)
    if(missing(filterId)){
        filterId <- deparse(subs)
        if(length(filterId)>1)
            filterId <- paste(gsub("^ *", "", filterId[2]), "...", sep="")
    }
    new("expressionFilter", filterId=filterId, expr=as.expression(subs),
        args=list(...), deparse=deparse(subs))
}

## construct expression by parsing a character string
char2ExpressionFilter <- function(expr, ..., filterId)
{
    subs <- parse(text=expr)
    if(missing(filterId))
         filterId <- expr
    new("expressionFilter", filterId=filterId, expr=subs,
        args=list(...), deparse=expr)
}



## ===========================================================================
## timeFilter
## ---------------------------------------------------------------------------
## Detect turbulences and abnormalities in the aquisition of flow data over
## time and gate them out. Argument 'bandwidth' sets the sensitivity, i.e.,
## the amount of local variance of the signal we want to allow. 'binSize'
## controls the size of the bins for the local variance and location
## estimation, 'timeParameter' can be used to explicitely give the paramter
## name of the time parameter (we will make an educated guess if this is not
## given)
## ---------------------------------------------------------------------------
setClass("timeFilter",
         representation(bandwidth="numeric",
                        binSize="numeric",
                        timeParameter="character"),
         contains="parameterFilter")

## contructor
timeFilter <- function(..., filterId="time", bandwidth=0.75,
                       binSize=NULL, timeParameter=NULL)
{
  x <- ..1
  if(is.list(x))
    pars <- unlist(x)
  else if (is.character(x) && length(x)>1)
    pars <- x
  else
    pars <- unlist(list(...))
  if(!all(is.character(pars)))
    stop("Only know how to deal with character data for the parameter ",
         "definition.", call.=FALSE)
  new("timeFilter", parameters=pars,
      bandwidth=bandwidth, binSize=as.numeric(binSize),
      timeParameter=as.character(timeParameter), filterId=filterId)
}




## ===========================================================================
## filterSet
## ---------------------------------------------------------------------------
## stores a list of filters from a gating sequence as an environment
## ---------------------------------------------------------------------------
setClass("filterSet",
         representation(env="environment"),
         prototype=prototype(env=new.env(hash=TRUE, parent=emptyenv())))

##constructor
filterSet <- function(...) {
    filters <- list(...)
    ## Allow the list(x, y, z) format as well.
    if(length(filters)==1 && is.list(filters[[1]]))
        filters <- filters[[1]]
    if(length(filters) == 0)
        new("filterSet", env=new.env(parent=emptyenv()))
    else
        as(filters, "filterSet")
}



## ===========================================================================
## filterReference
## ---------------------------------------------------------------------------
## References a filter (contained within a filterSet)
## ---------------------------------------------------------------------------
setClass("filterReference",
         representation(name="character", env="environment"),
         contains="filter")



## ===========================================================================
## setOperationFilter
## ---------------------------------------------------------------------------
setClass("setOperationFilter",
         representation(filters="list"),
         contains="concreteFilter")



## ===========================================================================
## unionFilter 
## ---------------------------------------------------------------------------
setClass("unionFilter",
         representation("setOperationFilter"))



## ===========================================================================
## intersectFilter 
## ---------------------------------------------------------------------------
setClass("intersectFilter",
         representation("setOperationFilter"))



## ===========================================================================
## complementFilter 
## ---------------------------------------------------------------------------
setClass("complementFilter",
         representation("setOperationFilter"),
         validity=function(object)
     { 
         if(length(object@filters) != 1) {
             warning("Complement filters can only operate on a ",
                     "single filter")
             return(FALSE)
         }
         TRUE
     })



## ===========================================================================
## subsetFilter 
## ---------------------------------------------------------------------------
setClass("subsetFilter",
         representation("setOperationFilter"),
         validity=function(object)
     {
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
## overlapping sets. The subset indices are stored as a matrix, where
## each row contains the results of a single filtering operation.
## ---------------------------------------------------------------------------
setClass("manyFilterResult",
         representation(subSet="matrix", dependency="ANY"),
         contains="filterResult")

##constructor
manyFilterResult <- function(filters, frameId, dependency=NULL)
{
    q <- new("manyFilterResult",
             filterDetails=sapply(filters, slot, "filterDetails"),
             subSet=do.call("cbind", lapply(filters, as, "logical")),
             dependency=dependency)
    colnames(q@subSet) <- sapply(filters, slot, "filterId")
    q
}



## ===========================================================================
## randomFilterResult
## ---------------------------------------------------------------------------
setClass("randomFilterResult",
         representation(subSet="numeric"),
         contains="filterResult")



## ===========================================================================
## filterSummary
## ---------------------------------------------------------------------------
setClass("filterSummary",
         representation(name="character", true="numeric",
                        count="numeric", p="numeric"))



## ===========================================================================
## transform
## ---------------------------------------------------------------------------
## Parameterize transforms so that we can describe them.
## ---------------------------------------------------------------------------
setClass("transform",
         representation(transformationId="character", "function"))

## Linear transform constructor
linearTransform <- function(transformationId, a=1, b=0)
{
    if(!is.double(a)) 
        stop("a must be numeric")
    if(!is.double(b))
        stop("b must be numeric")
    t <- new("transform", .Data=function(x)  x <- a*x+b)
    t@transformationId <- transformationId
    t
}

## Quadratic transformation constructor
quadraticTransform <- function(transformationId, a=1, b=1, c=0)
{
    if(!is.double(a)) 
        stop("a must be numeric")
    if(!is.double(b))
        stop("b must be numeric")
    if(!is.double(c))
        stop("c must be numeric")
    t <- new("transform", .Data=function(x) x <- a*x^2 + b*x + c)
    t@transformationId <- transformationId
    t
}

## Natural logarithm transformation constructor
lnTransform <- function(transformationId, r=1, d=1)
{
    if(!is.double(r) || r <= 0)
        stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
        stop("d must be numeric")
    t <- new("transform", .Data=function(x)
             x<-log(x)*(r/d))
    t@transformationId <- transformationId
    t
}

## Logarithm transformation constructor
logTransform <- function(transformationId, logbase=10, r=1, d=1)
{
    if(!is.double(r) || r <= 0)
        stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
        stop("d must be numeric")
    if(!is.double(r) || r <=0)
        stop("r must be numeric and positive")
    if(!is.double(logbase) || logbase <= 1)
        stop("logabse must be a pnumeric greater than 1")
    t <- new("transform", .Data=function(x) x <- log(x, logbase)*(r/d))
    t@transformationId <- transformationId
    t
}


## General biexponential transformation constructor
biexponentialTransform <-
    function(transformationId, a=.5, b=1, c=.5, d=1, f=0, w=0,
             tol=.Machine$double.eps^0.25, maxit=as.integer(5000))
{
    t <- new("transform", .Data=function(x)
             x <- .Call(biexponential_transform, x, a, b, c,
                        d, f, w, tol, maxit))
    t@transformationId <- transformationId
    t
}

## Logicle transformation constructor
logicleTransform <- function(transformationId, w=0, r=262144, d=5, ...)
{
    if(w>d)
        stop("Negative range decades must be smaller than total ",
             "number of decades")
    w <- w*log(10)
    d <- d*log(10)
    p <- if(w==0) 1 else uniroot(function(p) -w+2*p*log(p)/(p+1),
            c(.Machine$double.eps, 2*(w+d)))$root
    t <- new("transform", .Data=biexponentialTransform(transformationId,
                          a=r*exp(-(d-w)), b=1, c=r*exp(-(d-w))*p^2, d=1/p,
                          f=p^2-1, w=w, ...))
    t@transformationId <- transformationId
    t
}

## Truncation transformation constructor
truncateTransform <- function(transformationId, a=1)
{
    t <- new("transform", .Data=function(x){
        x[x<=a] <- a
        x
    })
    t@transformationId <- transformationId
    t
}

## Scale transformation constructor
scaleTransform <- function(transformationId, a=1, b=10^4)
{
    t <- new("transform", .Data=function(x) (x-a)/(b-a))
    t@transformationId <- transformationId
    t
}

## Split-scale transformation constructor
splitScaleTransform <- function(transformationId, maxValue=1023,
                                transitionChannel=64, r=192)
{
    maxChannel <- r + transitionChannel
    b <- transitionChannel/2
    d <- 2*log10(exp(1))*r/transitionChannel
    logt <- -2*log10(exp(1))*r/transitionChannel + log10(maxValue)
    t <- 10^logt
    a <- transitionChannel/(2*t)
    logCT <- (a*t+b)*d/r
    c <- 10^logCT/t
    tr <- new("transform", .Data= function(x){
        idx <- which(x <= t)
        idx2 <- which(x > t)
        if(length(idx2)>0)
            x[idx2] <- log10(c*x[idx2])*r/d
        if(length(idx)>0)
            x[idx] <- a*x[idx]+b
        x
    })
    tr@transformationId <- transformationId
    tr
}

## Hyperbolic Arcsin transformation constructor
arcsinhTransform <- function(transformationId, a=1, b=1, c=0)
{
    t <- new("transform", .Data=function(x) asinh(a+b*x)+c)
    t@transformationId <- transformationId
    t
}



## ===========================================================================
## parameterTransform
## ---------------------------------------------------------------------------
## A class that only applied a parameter transform to a subset of the
## parameters during an %on% operation.
## ---------------------------------------------------------------------------
setClass("parameterTransform",
         representation(parameters="character"),
         contains="transform")

## constructor
parameterTransform <- function(FUN, params)
    new("parameterTransform", .Data=as.function(FUN),
        parameters=as.character(params))



## ===========================================================================
## transformMap
## ---------------------------------------------------------------------------
## We want to be able to include transforms within a filter. First we need to
## know which parameters should be input filters
## ---------------------------------------------------------------------------
setClass("transformMap",
         representation(output="character", input="character", f="function"))



## ===========================================================================
## transformList
## ---------------------------------------------------------------------------
## A list of transformMaps
## ---------------------------------------------------------------------------
setClass("transformList",
         representation(transforms="list"),
         validity=function(object)
         if(all(sapply(object@transforms, is, "transformMap"))) TRUE else
         stop("All list items of a 'transformList' must be of class ",
              "'transformMap.'", call.=FALSE))

## constructor
transformList <- function(from, tfun, to=from)
{
    from <- unique(from)
    to <- unique(to)
    if(!is.character(from) || !is.character(to) || length(from) != length(to))
        stop("'from' and 'to' must be character vectors of equal length.",
             call.=FALSE)
    if(is.character(tfun))
        tfun <- lapply(tfun, get)
    if(!is.list(tfun)) tfun <- list(tfun)
    if(!all(sapply(tfun, is, "function")))
        stop("'tfun' must be a list of functions or a character vector ",
             "with the function names.", call.=FALSE)
    tfun <- rep(tfun, length(from))
    tlist <- mapply(function(x, y, z)
                    new("transformMap", input=x, output=y, f=z),
                    from, to, tfun[1:length(from)])
    return(as(tlist, "transformList"))
}



## ===========================================================================
## transformFilter
## ---------------------------------------------------------------------------
setClass("transformFilter",
         representation(transforms="transformList", filter="filter"),
         contains="concreteFilter")
