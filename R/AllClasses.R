## =========================================================================##
## =========================================================================##
##                    Class definitions and contructors                     ##
## =========================================================================##
## =========================================================================##






## ===========================================================================
##  Some helpers
## ---------------------------------------------------------------------------
## Check for the class of object x and its length and cast error if wrong
checkClass <- function(x, class, length=NULL, verbose=FALSE,
                       mandatory=TRUE)
{
    if(mandatory && missing(x))
        stop("Argument '", substitute(x), "' missing with no default",
             call.=verbose)
    msg <- paste("'", substitute(x), "' must be object of class '",
                 class, "'", sep="")
    fail <- !is(x, class)
    if(!is.null(length) && length(x) != length){
        fail <- TRUE
        msg <- paste(msg, "of length", length)
    }
    if(fail) stop(msg, call.=verbose) else invisible(NULL)     
}



## ===========================================================================
##  flowFrame
## ---------------------------------------------------------------------------
## A container for flow cytometry measurements with slots exprs, parameters
## and description. exprs contains measurement values, description contains 
## information from file headers of FCS file and parameters contains
## information about the FCS measurement parameters (i.e. channels) available.
## ---------------------------------------------------------------------------
setClass("flowFrame",                
         representation=representation(exprs="matrix",
         parameters="AnnotatedDataFrame",
         description="list"),
         prototype=list(exprs=matrix(numeric(0),
                        nrow=0,
                        ncol=0),
         parameters=new("AnnotatedDataFrame"),
         description=list(note="empty")))

## helper function to create empty AnnotatedDataFrame for the parameters slot
parDefault <- function(exp)
{
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
isValidParameters <- function(parameters, exprs)
{
    checkClass(parameters, "AnnotatedDataFrame")
    if(!all(c("name", "desc", "range", "minRange", "maxRange")
            %in% varLabels(parameters)))
        stop("The following columns are mandatory:\n  'name', 'desc',",
             "'range', 'minRange', 'maxRange'", call.=FALSE)
    if(!missing(exprs))
        if(!all(colnames(exprs) %in% parameters$name))
            stop("parameter description doesn't match colnames of the ",
                 "data matrix", call.=FALSE)
    return(TRUE)
}

## constructor
flowFrame <- function(exprs, parameters, description=list())
{
    if(!is.matrix(exprs) || !is.numeric(exprs) || is.null(colnames(exprs)))
        stop("Argument 'exprs' must be numeric matrix with colnames ",
             "attribute set", call.=FALSE)
    if(missing(parameters))
        parameters <- parDefault(exprs)
    else
        isValidParameters(parameters, exprs)
    checkClass(description, "list")
    new("flowFrame", exprs=exprs, parameters=parameters,
        description=description)
}



## ===========================================================================
##  flowSet
## ---------------------------------------------------------------------------
## A collection of several cytoFrames making up one experiment. Slots 
## frames, phenoData, colnames. Frames contains the cytoFrame objects,
## phenoData the experiment meta data and colnames the channel names.
## An additional character scalar _.name._ is stored in the environment
## which holds a name of the object that will be used in the workFlows.
## By storing it in the environment we don't have to add an additional
## slot and defunct old serialized flowSet objects.
## ---------------------------------------------------------------------------
setClass("flowSet",                   
         representation=representation(frames="environment",
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
flowSet <- function(..., phenoData, name)
{
    x <- list(...)
    if(length(x) == 1 && is.list(x[[1]]))
        x <- x[[1]]
    if(!all(sapply(x, is, "flowFrame")))
        stop("All additional arguments must be flowFrames")
    f <- as(x, "flowSet")
    if(!missing(phenoData))
        phenoData(f) <- phenoData
    if(!missing(name))
        identifier(f) <- name
    f
}



## ===========================================================================
## transform parent class and parameters
## ---------------------------------------------------------------------------
## Parameterize transforms so that we can describe them.
## ---------------------------------------------------------------------------
setClass("transform",
         representation=representation(transformationId="character",
                                       .Data="function"),
         prototype=prototype(transformationId=""))

setClass("parameters", contains="list")

setClassUnion("transformation", "transform")

setClassUnion("characterOrTransformation", c("character","transformation"))

setClassUnion("characterOrParameters", c("character","parameters"))

setClass("singleParameterTransform",
         representation=representation(parameters="transformation"),
         contains="transform")

setClass("nullParameter",
         representation=representation(dummy="numeric"))



## ===========================================================================
## Virtual filter and derived concreteFilter and parameterFilter
## ---------------------------------------------------------------------------
## A class describing a selection applied to a flow data matrix. Consist of
## a filterId and the names of the parameters to operate on (for parameter
## filters only). More specific filters all inherit from either of these two
## classes.
## ---------------------------------------------------------------------------
setClass("filter", 
         representation=representation("VIRTUAL",
         filterId="character"),
         prototype=prototype(filterId=""))

setClass("concreteFilter",
         contains="filter")

                                        # setClass("parameterFilter",
                                        #          representation=representation(parameters="character"),
                                        #          contains="concreteFilter",
                                        #          prototype=prototype(parameters=""))

setClass("parameterFilter", 
         representation(parameters="parameters"),
         contains="concreteFilter",
         prototype=prototype(parameters=new("parameters",.Data="NULL"))
         )


## ===========================================================================
## Rectangular gate
## ---------------------------------------------------------------------------
## A class describing a 2D rectangular region in the parameter space. Slots
## min and max hold the boundaries in the two dimensions.
## ---------------------------------------------------------------------------
setClass("rectangleGate",
         representation=representation(min="numeric",
         max="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultRectangleGate",
         min=-Inf,
         max=Inf)
         )

## parse '...' argument of a gate constructor. The return value is a list
## with parameters (as transforms) and values.
parseDots <- function(dl, collapseFirst=TRUE, len=NULL){
    parseItem <- function(i, x, len){
        ## We can return transforms directly
        y <- x[[i]]
        if(is(y, "transform")){
            dl[[i]] <<- NA
            y
        }else{
            li <- length(y)
            if(!is.character(y) && !is.null(len) && li!=len)
                stop("All additional arguments must be of length ",
                     len, call.=FALSE)
            if(!is.character(y) && li!=allLen)
                stop("All additional arguments must be of equal length ",
                     call.=FALSE)     
            if(!is.character(y) && is.null(names(x)[i]))
                stop("Additional arguments have to be named.",
                     call.=FALSE)
            if(is.character(y)){
                ## We return character scalars as unitytransforms
                dl[[i]] <<- NA
                unitytransform(y) 
            }else{
                ## For eerything else we make unitytransforms from the
                ## argument names
                unitytransform(names(x)[i])
            }
        }
    }
    ## We only parse ..1 if it is a list and drop all other arguments
    if(collapseFirst && length(dl) && is.list(dl[[1]]))
        dl <- dl[[1]]
    if(length(dl)){
        ## If ..1 is a character vector we return unitytransforms only
        if(is.character(dl[[1]]) && length(dl[[1]])>1)
            return(list(parameters=sapply(dl[[1]], unitytransform,
                        simplify=FALSE),
                        values=as.list(rep(NA, length(dl[[1]])))))
        ## If ..1 is a matrix we return unitytransforms and the matrix
        if(is.matrix(dl[[1]])){
            if(is.null(colnames(dl[[1]])))
                stop("Matrix of gate boundaries must have colnames.",
                     call.=FALSE)
            return(list(parameters=sapply(colnames(dl[[1]]), unitytransform,
                        simplify=FALSE), values=dl[[1]]))
        }
        ## All items in dl must be of equal length
        allLen <- if(is.character(dl[[1]])) length(dl[[min(length(dl), 2)]]) 
        else length(dl[[1]])
       
    }
    parms <- sapply(seq_along(dl), parseItem, dl, len, simplify=FALSE)
    return(list(parameters=parms, values=dl))
}

## Further process the output of parseDots to collapse individual arguments
prepareInputs <- function(parsed, .gate, ...)
{
    parms <- parsed$parameters
    values <- parsed$values
    if(missing(.gate)){
        if(any(sapply(values, is.na)))
            stop("The gate boundaries has to be provides as argument",
                 " '.gate'", call.=FALSE)
        if(!is.matrix(values)){
            sel <- sapply(values, is, "numeric")
            if(any(sel)){
                values <- matrix(sapply(values[sel], function(x){
                    if(length(x) ==2)
                        x <- sort(x)
                    x}), ncol=length(parms))
                parms <- parms[sel]
                colnames(values) <- sapply(parms, parameters)
            }
            return(list(parameters=parms, values=values))
        }
        if(!length(parms))
            stop("No arguments provided.", call.=FALSE)
        return(parsed)
    }else{
        if(is.matrix(.gate) && !is.null(colnames(.gate)))
            return(parseDots(list(.gate), ...))
        if(any(sapply(values, is.na))){
            if(ncol(.gate) != length(parms))
                stop("Number of parameters and dimensions of supplied",
                     " gate boundaries don't match.", call.=FALSE)
            return(list(parameters=parms, values=.gate))
        }
        if(!length(parms) || !all(sapply(parms, is, "unityTranform"))){
            return(prepareInputs(parseDots(list(.gate), ...)))
        }else{
            return(parsed)
        }
    }
}


## Constructor. We allow for the following inputs:
##  ... are named numerics, each of length 2
##  ... are transforms or a mix of transforms and characters, .gate is
##      the associated matrix of min and max values
##  ..1 is a named list of numerics
##  ..1 is a list of transformations or characters and .gate is the
##      associated matrix of min and max values, each of length 2
##  .gate is a matrix of min and max values with colnames = parameters
##  .gate is a named list of numerics, each of lenght 2
rectangleGate <- function(..., .gate, filterId="defaultRectangleGate")
{
    checkClass(filterId, "character", 1)
    parms <- parseDots(list(...), len=2)
    parms <- prepareInputs(parms, .gate, len=2)
    parms$values <- apply(parms$values, 2, sort)
    new("rectangleGate", filterId = filterId, parameters=parms$parameters,
        min=parms$value[1, ], max=parms$value[2, ])
}



## ===========================================================================
## Quadrant gate
## ---------------------------------------------------------------------------
## A class describing a gate which separates a 2D parameter space into
## four quadrants. Slot boundary holds a vector of length two indicating
## the quadrant boundaries in each of the two dimensions.
## ---------------------------------------------------------------------------
setClass("quadGate",
         representation=representation(boundary="numeric"),        
         contains="parameterFilter",
         prototype=list(filterId="defaultQuadGate",
         boundary=c(Inf, Inf)))

## Constructor. We allow for the following inputs:
##  ..1 and ..2 are named numerics of length 1
##  ..1 and ..2 are transforms or a mix of transforms and characters, .gate
##      is the associated numeric vector of boundary values of length 2
##  ..1 is a named list of numerics of length 1
##  ..1 is a list of transformations or characters and .gate is the
##      associated numeric vector of boundary values of length 2
##  .gate is a named list of numerics, each of lenght 1
quadGate <- function(..., .gate, filterId="defaultQuadGate")
{
    checkClass(filterId, "character", 1)
    if(!missing(.gate) && !is.list(.gate) && !is.matrix(.gate))
        .gate <- matrix(.gate, nrow=1)
    parms <- prepareInputs(parseDots(list(...), len=1), .gate, len=1)
    if(length(parms$parameters) !=2 || nrow(parms$value)!=1)
        stop("Expecting two named arguments or a single named vector\n",
             "of length 2 as input for gate boundaries.", call.=FALSE)
    new("quadGate", filterId=filterId, parameters=parms$parameters,
        boundary=as.numeric(parms$value))
}



## ===========================================================================
## Polygon gate
## ---------------------------------------------------------------------------
## A class describing a 2D polygonal region in the parameter space. Slot
## boundary holds the vertices of the polygon in a 2 colum matrix.
## ---------------------------------------------------------------------------
setClass("polygonGate",
         representation(boundaries="matrix"),
         contains="parameterFilter",
         prototype=list(filterId="defaultPolygonGate",
         boundaries=matrix(ncol=2, nrow=3)),
         validity=function(object)
     {
         msg <- TRUE
         if(!is.matrix(object@boundaries) || nrow(object@boundaries)<3 ||
            ncol(object@boundaries)!=2
            )
             msg <- paste("\nslot 'boundaries' must be a numeric matrix",
                          "of at least 3 rows and exactly 2 columns")
         return(msg)
     })

## Constructor. We allow for the following inputs:
##  ..1 and ..2 are named numerics, each of the same length
##  ..1 and ..2  are transforms or a mix of transforms and characters, .gate is
##      the associated matrix of polygon vertices of ncol=2
##  ..1 is a named list of numerics, each of the same length
##  ..1 is a list of transformations or characters and .gate is the
##      associated matrix of polygon vertices of ncol=2
##  .gate is a matrix of polygon vertices of ncol=2, colnames = parameters
##  .gate is a named list of two numerics, both of the same length
polygonGate <- function(..., .gate, boundaries, filterId="defaultPolygonGate")
{
    checkClass(filterId, "character", 1)
    if(missing(.gate))
        if(!missing(boundaries)){
            .Deprecated(msg=paste("The 'boundaries' argument is deprecated,",
                        "please use '.gate' instead."))
            .gate=boundaries
        }     
    parms <- prepareInputs(parseDots(list(...)), .gate)
    if(length(parms$parameters) !=2)
        stop("Polygon gates are only defined in two dimensions.",
             call.=FALSE)
    new("polygonGate", filterId=filterId, parameters=parms$parameters,
        boundaries=parms$values)
}



## ===========================================================================
## Polytope gate
## ---------------------------------------------------------------------------
## A class describing a nD polytope region in the parameter space. Slot a
## holds the coefficients of the linear equations for m halfspaces in n
## dimensions and b is a vector of m intercepts.
## ---------------------------------------------------------------------------
setClass("polytopeGate",
         representation(a="matrix",b="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultPolytopeGate", a=matrix(), b=1))

## Constructor. We allow for the following inputs:
##  b is always a numeric of length = ncol(a)
##  ... are named numerics, each of the same length
##  ...  are transforms or a mix of transforms and characters, a is
##      the associated matrix of coefficients
##  ..1 is a named list of numerics, each of the same length
##  ..1 is a list of transformations or characters and a is the
##      associated matrix of coefficients
##  .gate is a matrix of coefficients , colnames = parameters
##  .gate is a named list of numerics, all of the same length
polytopeGate <- function(..., .gate, b, filterId="defaultPolytopeGate")
{
    checkClass(filterId, "character", 1)
    checkClass(b, "numeric")
    parms <- prepareInputs(parseDots(list(...)), .gate)
    colnames(parms$values) <- sapply(parms$parameters, parameters)
    new("polytopeGate", filterId=filterId, parameters=parms$parameters,
        a=parms$values, b=b)
}



## ===========================================================================
## Ellipsoid gate
## ---------------------------------------------------------------------------
## A class describing an ellipsoid region in the parameter space. Slots
## mean and cov contain the mean values and the covariance matrix describing
## the ellipse, slot distance holds a scaling factor, i.e., the Mahalanobis
## distance.
## ---------------------------------------------------------------------------
setClass("ellipsoidGate",
         representation(mean="numeric",
                        cov="matrix",
			distance="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultEllipsoidGate",
         mean=numeric(), cov=matrix(), distance=1),
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

## Constructor. We allow for the following inputs:
##  mean always is a numeric of the same length as number of dimensions,
##  distance is always a vector of length 1
##  ... are named numerics, each of the same length
##  ...  are transforms or a mix of transforms and characters, .gate is
##      the associated covariance matrix
##  ..1 is a named list of numerics, each of the same length
##  ..1 is a list of transformations or characters and .gate is the
##      associated covariance matrix
##  .gate is the covariance matrix, colnames=parameters
##  .gate is a named list of numerics, each of the same length
ellipsoidGate <- function(..., .gate, mean, distance=1,
                          filterId="defaultEllipsoidGate")
{
    checkClass(filterId, "character", 1)
    checkClass(mean, "numeric")
    checkClass(distance, "numeric", 1)
    parms <- prepareInputs(parseDots(list(...)), .gate)
    names(mean) <- sapply(parms$parameters, parameters)
    new("ellipsoidGate", filterId=filterId, parameters=parms$parameters,
        cov=parms$values, mean=mean, distance=distance)
}



## ===========================================================================
## norm2Filter
## ---------------------------------------------------------------------------
## A class to describe the fit of a bivariate normal distribution.
## Slot method is a character describing the method used to compute the
## covariance matrix, slot scale.factor holds a numeric representing the
## Mahalanobis distance. Slot transformation holds a list of length
## giving transformations, if applicable that are applied to the data
## before gating. n is the number of points used in the subsampling step.
## ---------------------------------------------------------------------------
setClass("norm2Filter",
         representation=representation(method="character",
         scale.factor="numeric",
         n="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultNorm2Filter",
         scale.factor=1,
         transformation=list(),
         method="covMcd",
         n=50000))

## Constructor. We allow for the following inputs:
##  method is always a character and scale.factor and n both are always
##     numerics, all of length 1
##  x and y are characters of length 1 or a mix of characters and
##     transformations
##  x is a character of length 2 and y is missing
##  x is a list of characters and/or transformations, y is missing
norm2Filter <- function(x, y, method="covMcd", scale.factor=1,
                        n=50000, filterId="defaultNorm2Filter")
{
    checkClass(method, "character", 1)
    checkClass(scale.factor, "numeric", 1)
    checkClass(n, "numeric", 1)
    checkClass(filterId, "character", 1)
    if(missing(y)) {
        if(length(x)==1)
            stop("You must specify two parameters for a norm2 gate.")
        if(length(x)>2)
            warning("Only the first two parameters will be used.")
        y=x[[2]]
        x=x[[1]]
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
         representation=representation(populations="character"),
         prototype=list(filterId="defaultKmeansFilter"),
         contains="parameterFilter")

## Constructor. We allow for the following inputs:
##  ..1 is transform and .2 is some vector that can be coerced to character
##  ..1 is some vector that can be coerced to character
kmeansFilter <- function(..., filterId="defaultKmeansFilter")
{
    checkClass(filterId, "character", 1)
    ll <- list(...)
    if(length(ll)){
        n <- names(ll)[1]
        if(is.list(ll[[1]])){
            n <- names(ll[[1]])[1]
            ll[[1]] <- unlist(ll[[1]], recursive=FALSE)
        }
        parameter <- if(is(ll[[1]], "transform")) ll[[1]] else n
        populations <- if(is(parameter, "transform")){
            if(length(ll)==1)
                stop("List of populations needs to be provided as ",
                     "an additional argument.", call.=FALSE) 
            as.character(unlist(ll[[2]]))} else as.character(unlist(ll[[1]]))
    }else{
        stop("No arguments provided.", .call=FALSE)
    }
    new("kmeansFilter", parameters=parameter,
        populations=populations, filterId=filterId)
}



## ===========================================================================
## curv1Filter
## ---------------------------------------------------------------------------
## This filter can hold parameters to find siginficant high density regions
## in one dimension based on Matt Wand's feature software. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("curv1Filter",
         representation=representation(bwFac="numeric",
         gridsize="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultCurv1Filter",
         bwFac=1.2,
         gridsize=rep(151, 2)))

## Constructor. We allow for the following inputs:
##  bwFac is always a numeric of length 1 and gridsize is always a numeric
##     of length 2
##  x is either a character or a transformation
curv1Filter <- function(x, bwFac=1.2, gridsize=rep(151, 2),
                        filterId="defaultCurv1Filter")
{
    checkClass(filterId, "character", 1)
    checkClass(bwFac, "numeric", 1)
    checkClass(gridsize, "numeric", 2)
    new("curv1Filter", parameters=x, bwFac=bwFac,
        gridsize=gridsize, filterId=as.character(filterId))
}



## ===========================================================================
## curv2Filter
## ---------------------------------------------------------------------------
## This filter can hold parameters to find siginficant high density regions
## in two dimensions based on Matt Wand's feature software. This generates a
## multipleFilterResult
## ---------------------------------------------------------------------------
setClass("curv2Filter",
         representation=representation(bwFac="numeric",
         gridsize="numeric"),
         contains="parameterFilter",
         prototype=list(filterId="defaultCurv2Filter",
         bwFac=1.2,
         gridsize=rep(151, 2)))

## Constructor. We allow for the following inputs:
##  bwFac is always a numeric of length 1 and gridsize is always a numeric
##     of length 2
##  x and y are characters of length 1 or a mix of characters and
##     transformations
##  x is a character of length 2 and y is missing
##  x is a list of characters and/or transformations, y is missing
curv2Filter <- function(x, y, filterId="defaultCurv2Filter",
                        bwFac=1.2, gridsize=rep(151, 2))
{
    checkClass(filterId, "character", 1)
    checkClass(bwFac, "numeric", 1)
    checkClass(gridsize, "numeric", 2) 
    if(missing(y)) {
        if(length(x)==1)
            stop("You must specify two parameters for a curv2Filter.",
                 call.=FALSE)
        if(length(x)>2)
            warning("Only using parameters '", x[1], "' and '", x[2],
                    "'.", call.=FALSE)
        y=x[[2]]
        x=x[[1]]
    }
    new("curv2Filter", parameters=list(x,y), bwFac=bwFac,
        gridsize=gridsize, filterId=as.character(filterId))
}



## ===========================================================================
## sampleFilter
## ---------------------------------------------------------------------------
## Sample 'size' rows from a flowFrame. 
## ---------------------------------------------------------------------------
setClass("sampleFilter",
         representation=representation(size="numeric"),
         contains="concreteFilter",
         prototype=list(size=10000, filterId="defaultSampleFilter"))

##Constructor: We allow for the following inputs:
##  size is always a numeric of length 1
sampleFilter <- function(size, filterId="defaultSampleFilter")
{
    checkClass(filterId, "character", 1)
    checkClass(size, "numeric", 1)
    new("sampleFilter", filterId=filterId, size=size)
}



## ===========================================================================
## expressionFilter
## ---------------------------------------------------------------------------
## Let's us encapsulate an expression as a gate. There also is a constructor
## to create the filter from a character representation of the expression
## which is helpful for programmatic use. The args slot can contain additional
## arguments that are passed on to the evaluation environment. deparse stores
## a deparsed version of the expression.
## ---------------------------------------------------------------------------
setClass("expressionFilter",
         representation=representation(expr="expression",
         args="list",
         deparse="character"),
         contains="concreteFilter",
         prototype=list(filterId="defaultExpressionFilter",
         exprs=expression(rep(TRUE, length(get(ls()[1])))),
         args=list(),
         deparse="default"))

## Constructor: We allow for the following inputs:
##  expr is always an expression
##  ... are further arguments to the expression
expressionFilter <- function(expr, ..., filterId="defaultExpressionFilter")
{
    subs <- substitute(expr)
    if(missing(filterId)){
        filterId <- deparse(subs)
        if(length(filterId)>1)
            filterId <- paste(gsub("^ *", "", filterId[2]), "...", sep="")
    }else checkClass(filterId, "character", 1)
    new("expressionFilter", filterId=filterId, expr=as.expression(subs),
        args=list(...), deparse=deparse(subs))
}

## Constructor from a character string: We allow for the following inputs:
##  expr is always a character string
char2ExpressionFilter <- function(expr, ...,
                                  filterId="defaultExpressionFilter")
{
    checkClass(expr, "character", 1)
    subs <- parse(text=expr)
    if(missing(filterId))
        filterId <- expr
    else
        checkClass(filterId, "character", 1)
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
## given).
## ---------------------------------------------------------------------------
setClass("timeFilter",
         representation=representation(bandwidth="numeric",
         binSize="numeric",
         timeParameter="character"),
         contains="parameterFilter",
         prototype=list(filterId="defaultTimeFilter",
         bandwidth=0.75,
         binSize=NULL,
         timeParameter=NULL))

## Constructor: We allow for the following inputs:
##  bandwidth and binSize are always numerics of lenght 1, timeParameter
##      is always a character of length 1
##  ..1 is a character
##  ..1 is a list of character and/or transformations
##  ... are characters and/or transformations
timeFilter <- function(..., bandwidth=0.75, binSize, timeParameter,
                       filterId="defaultTimeFilter")
{
    checkClass(bandwidth, "numeric", 1)
    
    if(!missing(binSize))
        checkClass(binSize, "numeric", 1)
    else
        binSize <- NULL
    if(!missing(timeParameter))
        checkClass(timeParameter, "character", 1)
    else
        timeParameter <- NULL
    checkClass(filterId, "character", 1)
    parms <- parseDots(list(...))
    new("timeFilter", parameters=parms$parameters,
        bandwidth=bandwidth, binSize=as.numeric(binSize),
        timeParameter=as.character(timeParameter), filterId=filterId)
}



## ===========================================================================
## filterSet
## ---------------------------------------------------------------------------
## Stores a list of filters from a gating sequence as an environment. Although
## this will be kept as part of flowCore we encourage to use the new workflow
## infrastructure instead.
## ---------------------------------------------------------------------------
setClass("filterSet",
         representation=representation(env="environment",
         name="character"),
         prototype=prototype(env=new.env(hash=TRUE, parent=emptyenv()),
         name="Filter Set"))

## constructor
filterSet <- function(..., name="default") {
    filters <- list(...)
    ## Allow the list(x, y, z) format as well.
    if(length(filters)==1 && is.list(filters[[1]]))
        filters <- filters[[1]]
    if(length(filters) == 0)
        new("filterSet", env=new.env(parent=emptyenv()), name=name)
    else{
        tmp <- as(filters, "filterSet")
        tmp@name <- name
        tmp
    }
    
}



## ===========================================================================
## filterReference
## ---------------------------------------------------------------------------
## References a filter (contained within a filterSet). Everything is just
## passed to the referenced filter. This may be better handled by the type
## system by having "real" filters inherit from concreteFilter (or something)
## and then simply having a setAs(), but I think that will be too much work
## for filter authors.
## ---------------------------------------------------------------------------
setClass("filterReference",
         representation=representation(name="character",
         env="environment"),
         contains="filter")

## Constructor from an environment
setMethod("filterReference",
          signature("environment", "character"),
          function(from, name) {
              new("filterReference", name=name, env=from)
          })

## Constructor from another filterSet
setMethod("filterReference",
          signature("filterSet", "character"),
          function(from,name)
      {
          new("filterReference", env=from@env,
              name=name)
      })



## ===========================================================================
## setOperationFilter
## ---------------------------------------------------------------------------
## Superclass for union intersect, complement and subset filter, which all
## consist of two or more component filters
## ---------------------------------------------------------------------------
setClass("setOperationFilter",
         representation=representation(filters="list"),
         contains="concreteFilter")



## ===========================================================================
## unionFilter 
## ---------------------------------------------------------------------------
## The union of two filters, .i.e, the logical | operation.
## A simple optimization would be to linearize the union of a filter and
## another union filter.
## ---------------------------------------------------------------------------
setClass("unionFilter",
         representation=representation("setOperationFilter"))

## constructor from two filters
setMethod("|",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("unionFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "or", identifier(e2)))
      })

## constructor from a list of filters and a filter and vice versa
setMethod("|",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "|", e2=e2))
setMethod("|",
          signature=signature(e1="filter",
          e2="list"),
          definition=function(e1, e2) lapply(e2, "|", e1=e1))



## ===========================================================================
## intersectFilter 
## ---------------------------------------------------------------------------
## The intersection of two filters, i.e, the logical & operation.
## This is somewhat different from the %subset% operation because
## some filters depend on the data and would return different results
## when applied to the full dataset.
## --------------------------------------------------------------------------
setClass("intersectFilter",
         representation=representation("setOperationFilter"))

## constructor from two filters
setMethod("&",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("intersectFilter", filters=list(e1, e2),
              filterId=paste(identifier(e1), "and", identifier(e2)))
      })

## constructor from a list of filters and a filter and vice versa
setMethod("&",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "&", e2=e2))
setMethod("&",
          signature=signature(e1="filter",
          e2="list"),
          definition=function(e1, e2) lapply(e2, "&", e1=e1))



## ===========================================================================
## complementFilter 
## ---------------------------------------------------------------------------
## The complement of a filters, i.e, the logical ! operation.
## ---------------------------------------------------------------------------
setClass("complementFilter",
         representation=representation("setOperationFilter"),
         validity=function(object)
     { 
         if(length(object@filters) != 1) {
             warning("Complement filters can only operate on a ",
                     "single filter")
             return(FALSE)
         }
         TRUE
     })


## constructor
setMethod("!",
          signature=signature(x="filter"),
          definition=function(x)
      {
          new("complementFilter",filters=list(x),
              filterId=paste("not",identifier(x)))
      })



## ===========================================================================
## subsetFilter 
## ---------------------------------------------------------------------------
## Combining two filters in a way that the RHS filter  takes the subset
## of the LHS filter as input. For many cases this is equivalent to an
## intersection filter, the only differnce is in data-driven filters.
## ---------------------------------------------------------------------------
setClass("subsetFilter",
         representation=representation("setOperationFilter"),
         validity=function(object)
     {
         if(length(object@filters) != 2) {
             warning("Subset filters are only defined as binary operators")
             return(FALSE)
         }
         TRUE
     })

## constructor from two filters. %&% is an alias for %subset%
setMethod("%subset%",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2)
      {
          new("subsetFilter",
              filters=list(e1, e2), filterId=paste(identifier(e1),"in",
                                    identifier(e2)))
      })
setMethod("%&%",
          signature=signature(e1="filter",
          e2="filter"),
          definition=function(e1, e2) e1 %subset% e2)

## constructor from a list of filters and a filter
setMethod("%subset%",
          signature=signature(e1="list",
          e2="filter"),
          definition=function(e1, e2) lapply(e1, "%subset%", e2=e2))

## constructor from a filterSet and a filter
setMethod("%subset%",
          signature=signature(e1="filterSet",
          e2="filter"),
          definition=function(e1,e2)
      {
          ## Make a copy of the filterSet, preserving R semantics
          x <- as(as(e1, "list"), "filterSet")
          n <- names(e1)
          x[[""]] <- e2
          target <- as.symbol(identifier(e2))
          for(i in n){
              x[[""]] <- substitute(~ a %subset% b, list(a=as.symbol(i),
                                                         b=target))
          }
          x
      })



## ===========================================================================
## filterResult
## ---------------------------------------------------------------------------
## A container for the results after applying a filter to flow cytometry
## data with slots frameId (identifier of the object) and filterDetails,
## which is a list containing and further describing the input filter.
## ---------------------------------------------------------------------------
setClass("filterResult",
         representation=representation(frameId="character",
         filterDetails="list"),
         contains="concreteFilter",
         prototype=list(frameId="Filter Result",
         filterDetails=list()))



## ===========================================================================
## logicalFilterResult
## ---------------------------------------------------------------------------
## Resuls from a filtering operation that only produces a single population.
## Slot subet is a logical vector indicating the population membership of the
## data in the gated flowFrame.
## ---------------------------------------------------------------------------
setClass("logicalFilterResult",
         representation=representation(subSet="logical"),
         contains="filterResult")



## ===========================================================================
## multipleFilterResult
## ---------------------------------------------------------------------------
## Resuls from a filtering operation that produces multiple populations.
## Slot subet is a factor vector indicating the population membership of the
## data in the gated flowFrame. Factor names are used as population names.
## ---------------------------------------------------------------------------
setClass("multipleFilterResult",
         representation=representation(subSet="factor"),
         contains="filterResult")



## ===========================================================================
## manyFilterResult
## ---------------------------------------------------------------------------
## A special case of multipleFilterResult that arises when there are
## overlapping sets. The subset indices are stored as a matrix, where
## each row contains the results of a single filtering operation.
## ---------------------------------------------------------------------------
setClass("manyFilterResult",
         representation=representation(subSet="matrix",
         dependency="ANY"),
         contains="filterResult")

##constructor
manyFilterResult <- function(filters, frameId, dependency=NULL)
{
    q <- new("manyFilterResult",
             filterDetails=lapply(filters, slot, "filterDetails"),
             subSet=do.call("cbind", lapply(filters, as, "logical")),
             dependency=dependency)
    colnames(q@subSet) <- sapply(filters, slot, "filterId")
    q
}



## ===========================================================================
## randomFilterResult
## ---------------------------------------------------------------------------
## A result of a filtering operation where the population membership is
## considered to be stochastic rather than absolute. Currently there is no
## implementation of a filter that produces such a filterResult, although
## norm2Filter, curvFilters and the t-mixture filters in flowClust are
## obvious candidates.
## ---------------------------------------------------------------------------
setClass("randomFilterResult",
         representation=representation(subSet="numeric"),
         contains="filterResult")



## ===========================================================================
## filterResultList
## ---------------------------------------------------------------------------
## A list of filterResults which typically is generated when applying a
## filter to a whole flowSet. This is a class union of list and filterResult
## and mainly exists to allow for method dispatch and sanity checking.
## FIXME: Do we want to allow for mixed filter classes in the list?
## ---------------------------------------------------------------------------
setClass("filterResultList",
         contains=c("list", "filterResult"))

## Check if a filterResultList matches a flowSet. If strict=TRUE, the
## function will also check whether all items in the filterResultSet
## are of equal type and produce the same number of populations.
validFilterResultList <- function(fres, set, strict=TRUE)
{
    res <- TRUE
    checkClass(fres, "filterResultList")
    checkClass(strict, "logical", 1)
    if(!missing(set)){
        checkClass(set, "flowSet")
        if(res <- !all(names(fres) == sampleNames(set)))
            warning("Sample names don't match between flowSet and ",
                    "filterResultList", call.=FALSE)
    }
    if(strict){
        fTypes <- sapply(fres, function(x) class(x))
        if(length(unique(fTypes)) != 1){
            warning("Not all filterResults in the list are of equal",
                    " type.", call.=FALSE)
            res <- FALSE
        }
        nrPops <- sapply(fres, function(x) length(x))
        if(length(unique(nrPops)) != 1){
            warning("Not all filterResults in the list share the",
                    " same number of sub-populations.", call.=FALSE)
            res <- FALSE
        }
        return(res)
    }
}


## ===========================================================================
## filterSummary
## ---------------------------------------------------------------------------
## A class containing the results of calling summary methods on filterResult.
## In the case of multipleFilterrResults, the individual slots(except 'count')
## will be vectors.
## Slots are:
##   - name:  The name of the summary, usually this will be set to be the
##            identifier of the filterResult, or the names of the individual
##            populations for a multipleFilterResult
##   - true:  The number of events in the filter (or the individual
##            populations)
##   - count: The total number of events the filter was applied on
##   - p:     The ratio of events within the filter (i.e., true/count)
## ---------------------------------------------------------------------------
setClass("filterSummary",
         representation=representation(name="character",
         true="numeric",
         count="numeric",
         p="numeric"))



## ===========================================================================
## filterSummaryList
## ---------------------------------------------------------------------------
## A list of filterSummaries which typically is generated when summarizing a
## filterResultList. This directly extends the list class  and mainly exists
## to allow for method dispatch.
## ---------------------------------------------------------------------------
setClass("filterSummaryList",
         contains="list")



## ===========================================================================
## transform functions
## ---------------------------------------------------------------------------
## Constructors for the different varieties of transforms. All of these
## create objects of the basic class 'transform', unless stated otherwise.
## ---------------------------------------------------------------------------
## linear (polynomial) transform constructor
linearTransform <- function(transformationId, a=1, b=0)
{
    checkClass(a, "numeric")
    checkClass(b, "numeric")
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
## A class used to map parameters of a transform during %on% operations.
## ---------------------------------------------------------------------------
setClass("parameterTransform",
         representation=representation(parameters="character"),
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
         representation=representation(output="character",
         input="character",
         f="function"))



## ===========================================================================
## transformList
## ---------------------------------------------------------------------------
## A list of transformMaps
## ---------------------------------------------------------------------------
setClass("transformList",
         representation=representation(transforms="list"),
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
## FIXME: I have no clue what that is supposed to be but my guess is that it
## can go away once we have the new transformations in place
## ---------------------------------------------------------------------------
setClass("transformFilter",
         representation=representation(transforms="transformList",
         filter="filter"),
         contains="concreteFilter")



## ===========================================================================
## compensation
## ---------------------------------------------------------------------------
## A class to define a compensation operation.
## Slots are:
##   - compensationId: The identifier of the object
##   - spillover:      The spillover matrix
##   - parameters:     The parameters for which the data is to be compensated,
##                     an object of class parameters 
## ---------------------------------------------------------------------------
setClass("compensation",
         representation(spillover="matrix",
                        compensationId="character",
                        parameters="parameters"),
         prototype=prototype(spillover=matrix(),
                             compensationId="default",
                             parameters=new("parameters",.Data="")))

## Constructor: We allow for the following inputs:
##  spillover is always a symmetric numerical matrix with colnames set
##  invert is always a logical of length 1
##  ..1 is a character vector
##  ..1 is a list of character and/or transformations
##  ..1 is a matrix and spillover is missing
##  ... are characters and/or transformations
## If parameters are given explicitely they need to match the colnames
## of the spillover matrix.
compensation <- function(..., spillover, inv=FALSE,
                         compensationId="defaultCompensation")
{
    parms <- parseDots(list(...))
    if(missing(spillover))
        spillover <- as.matrix(parms$values)
    if(!is.matrix(spillover) || !is.numeric(spillover) ||
       ncol(spillover) != nrow(spillover))
        stop("'spillover' must be numeric matrix with same number of ",
             "rows and columns", call.=FALSE)
    if(is.null(colnames(spillover)))
        stop("Spillover matrix must have colnames", call.=FALSE)
    checkClass(compensationId, "character", 1)
    checkClass(inv, "logical", 1)
    if(!length(parms$parameters))
        parms <- sapply(colnames(spillover), unitytransform)

 #   if(all(sapply(is(parms$parameters, "unitytransform"))) && !all(sapply(parms$parameters, parameters) %in% colnames(spillover)))
 if(all(sapply(parms$parameters,function(x){ is(x,"unitytransform")}
              )
        ) && !all(sapply(parms$parameters, parameters) %in% colnames(spillover)))
        
    stop("Parameters and column names of the spillover matrix ",
             "don't match.", call.=FALSE)
    if(inv)
        spillover <- solve(spillover/max(spillover))
    new("compensation", spillover=spillover, 
        compensationId=compensationId,
        parameters=new("parameters", parms$parameters))
}



## ===========================================================================
## compensatedParameter
## ---------------------------------------------------------------------------
## FIXME NG: Please document
## ---------------------------------------------------------------------------
setClass("compensatedParameter",
          contains=c("transform"),
          representation=representation(parameters="character",spillRefId="character",
                                        searchEnv="environment"
                                       )
        )

## Constructor
compensatedParameter <- function(parameters,
                                 spillRefId="defaultCompensatedParameter",
                                 transformationId="defaultTransformationId",
                                 searchEnv)
{
    
    new("compensatedParameter", parameters=parameters, spillRefId=spillRefId,
        transformationId=transformationId,searchEnv=searchEnv)
}
    	 

## ===========================================================================
## fcReference
## ---------------------------------------------------------------------------
## The general definition of a reference: We have an identifier and an
## evaluation environment to do the lookup in. In the context of workflows,
## the environment will live inside the workFlow object.
## Rather than storing the actual items, the workflow framework will provide
## a more reference based semantic. To this end, we define a number
## of different reference classes to be used in the subsequent class
## definition. Objects of class "fcNullReference" allow to assign
## empty (hence unresolvable) references without breaking dispatch
## or class validity. The taks of creating the correct reference type is
## handled by the workFlow-specific assign methods, however these only call
## the appropriate fcReference constructors, so those need to make sure that
## all necesary side-effects take place.
## ---------------------------------------------------------------------------
## Create quasi-random guids. This is only based on the time stamp,
## not on MAC address or similar.
guid <- function()
    as.vector(format.hexmode(as.integer(Sys.time())/
                             runif(1)*proc.time()["elapsed"]))
setClass("fcReference", 
         representation=representation(ID="character",
         env="environment"),
         prototype=prototype(ID=paste("ref", guid(), sep="_"),
         env=new.env(parent=emptyenv()))
         )

## Set a value in the alias table of the workFlow 
setAlias <- function(alias, value, workflow)
{
    checkClass(alias, "character", 1)
    checkClass(value, "character", 1)
    checkClass(workflow, "workFlow")
    workflow <- alias(workflow)
    workflow[[alias]] <- unique(c(workflow[[alias]], value))
    return(invisible(NULL))
}

## Get a value from the alias table of the workFlow 
getAlias <- function(alias, workflow)
{
    checkClass(alias, "character")
    checkClass(workflow, "workFlow")
    fun <- function(x)
        if(x %in% ls(workflow)) x else alias(workflow)[[x]]
    return(as.vector(sapply(alias, fun)))
}

## remove alias for an identifier
rmAlias <- function(value, workflow)
{
    checkClass(value, "character", 1)
    checkClass(workflow, "workFlow")
    workflow <- alias(workflow)
    ind <- names(which(sapply(as.list(workflow), function(x)
                              value %in% x)==TRUE))
    for(i in ind){
        tmp <- workflow[[i]]
        tmp <- setdiff(tmp, value)
        if(!length(tmp))
            rm(list=i, envir=workflow)
        else
            workflow[[i]] <- tmp
    }
    return(invisible(NULL))
}

## Figure out which reference type to create, based on the class of 'value'
refType <- function(value)
{
    if(is(value, "flowFrame") || is(value, "flowSet")) "fcDataReference"
    else if(is(value, "filterResult")) "fcFilterResultReference"
    else if(is(value, "filter")) "fcFilterReference"
    else if(is(value, "actionItem")) "fcActionReference"
    else if(is(value, "view")) "fcViewReference"
    else if(is(value, "compensation")) "fcCompensateReference"
    else if(is(value, "transformList")) "fcTransformReference"
    else if(is(value, "graphNEL")) "fcTreeReference"
    else if(is(value, "environment")) "fcAliasReference"
    else if(is.null(value)) "fcNullReference"
    else "fcReference"
}

## Create useful identifiers for references
refName <- function(value)
{
    prefix <- if(is(value, "flowFrame") || is(value, "flowSet")) "dataRef"
    else if(is(value, "filterResult")) "fresRef"
    else if(is(value, "filter")) "filterRef"
    else if(is(value, "actionItem")) "actionRef"
    else if(is(value, "view")) "viewRef"
    else if(is(value, "compensation")) "compRef"
    else if(is(value, "transformList")) "transRef"
    else if(is(value, "graphNEL")) "treeRef"
    else if(is(value, "environment")) "aliasRef"
    else if(is.null(value)) "nullRef"
    else "genericRef"
    return(paste(prefix, guid(), sep="_"))
}

## constructor
fcReference <- function(ID=paste("genericRef", guid(), sep="_"),
                        env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "environment")
    ref <- new("fcReference", ID=ID, env=env)
    setAlias(substitute(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcStructureReference
## ---------------------------------------------------------------------------
## This only exists to subclass structural parts of the workFlow like views
## and actionItems to allow for group-wise method dispatch.
## ---------------------------------------------------------------------------
setClass("fcStructureReference",
         contains=list("VIRTUAL",
         "fcReference"))



## ===========================================================================
## fcTreeReference
## ---------------------------------------------------------------------------
## A reference to a graphNEL object representing the workflow tree
## ---------------------------------------------------------------------------
setClass("fcTreeReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("treeRef", guid(), sep="_"))
         )

## constructor
fcTreeReference <- function(ID=paste("treeRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcTreeReference", ID=ID, env=env@env)
    setAlias("tree", identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcAliasReference
## ---------------------------------------------------------------------------
## A reference to an environment object representing the alias table
## ---------------------------------------------------------------------------
setClass("fcAliasReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("aliasRef", guid(), sep="_"))
         )

## constructor
fcAliasReference <- function(ID=paste("aliasRef", guid(), sep="_"),
                             env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    new("fcAliasReference", ID=ID, env=env@env)
}



## ===========================================================================
## fcDataReference
## ---------------------------------------------------------------------------
## A reference to a data type object (a flowSet or flowFrame)
## ---------------------------------------------------------------------------
setClass("fcDataReference",
         contains="fcReference",
         prototype=prototype(ID=paste("dataRef", guid(), sep="_"))
         )

## constructor
fcDataReference <- function(ID=paste("dataRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcDataReference", ID=ID, env=env@env)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcActionReference
## ---------------------------------------------------------------------------
## A reference to an action item object (a gate, transformation or
## compensation)
## ---------------------------------------------------------------------------
setClass("fcActionReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("actionRef", guid(), sep="_"))
         )

## constructor
fcActionReference <- function(ID=paste("actionRef", guid(), sep="_"),
                              env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcActionReference", ID=ID, env=env@env)
    setAlias(names(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcViewReference
## ---------------------------------------------------------------------------
## A reference to a view object. This allows to bind action items to
## particular views on the data
## ---------------------------------------------------------------------------
setClass("fcViewReference",
         contains="fcStructureReference",
         prototype=prototype(ID=paste("viewRef", guid(), sep="_"))
         )

## constructor
fcViewReference <- function(ID=paste("viewRef", guid(), sep="_"),
                            env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcViewReference", ID=ID, env=env@env)
    setAlias(names(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcFilterResultReference
## ---------------------------------------------------------------------------
## A reference to a filterResult object. We need this to store filterResults
## along with the views the filtering generated, without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
setClass("fcFilterResultReference",
         contains="fcReference",
         prototype=prototype(ID=paste("fresRef", guid(), sep="_"))
         )

## constructor
fcFilterResultReference <- function(ID=paste("fresRef",
                                    guid(), sep="_"),
                                    env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcFilterResultReference", ID=ID, env=env@env)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcFilterReference
## ---------------------------------------------------------------------------
## A reference to a filter object. We need this to store filters within 
## a gateActionItem without unnecessarily copying things.
## ---------------------------------------------------------------------------
setClass("fcFilterReference",
         contains="fcReference",
         prototype=prototype(ID=paste("filterRef", guid(), sep="_"))
         )

## constructor
fcFilterReference <- function(ID=paste("filterRef",
                              guid(), sep="_"),
                              env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcFilterReference", ID=ID, env=env@env)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcCompensateReference
## ---------------------------------------------------------------------------
## A reference to a compensation object. We need this to store a compensation
## within a compensateActionItem without unnecessarily copying things.
## ---------------------------------------------------------------------------
setClass("fcCompensateReference",
         contains="fcReference",
         prototype=prototype(ID=paste("compRef", guid(), sep="_"))
         )

## constructor
fcCompensateReference <- function(ID=paste("compRef",
                                  guid(), sep="_"),
                                  env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcCompensateReference", ID=ID, env=env@env)
    setAlias(identifier(get(ref)), identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcTransformReference
## ---------------------------------------------------------------------------
## A reference to a transformation object. For now we use transformList
## objects until transformation is a more usefull class. We need this to
## store a trasnforamtion within a transformActionItem without unnecessarily
## copying things.
## ---------------------------------------------------------------------------
setClass("fcTransformReference",
         contains="fcReference",
         prototype=prototype(ID=paste("transRef", guid(), sep="_"))
         )

## constructor
fcTransformReference <- function(ID=paste("transRef",
                                 guid(), sep="_"),
                                 env=new.env(parent=emptyenv()))
{
    checkClass(ID, "character", 1)
    checkClass(env, "workFlow")
    ref <- new("fcTransformReference", ID=ID, env=env@env)
    setAlias("no scheme yet", identifier(ref), env)
    return(ref)
}



## ===========================================================================
## fcNullReference
## ---------------------------------------------------------------------------
## A NULL reference to be used whenever a slot is supposed to be empty
## ---------------------------------------------------------------------------
setClass("fcNullReference",
         contains=c("fcDataReference",
         "fcActionReference",
         "fcViewReference",
         "fcFilterResultReference",
         "fcFilterReference",
         "fcCompensateReference",
         "fcTransformReference",
         "fcTreeReference",
         "fcAliasReference"),
         prototype=prototype(ID=paste("nullRef", guid(), sep="_"))
         )

fcNullReference <- function(...) new("fcNullReference")



## ===========================================================================
## workFlow
## ---------------------------------------------------------------------------
## This class is intended to store all necessary information about a
## workflow. Since the individual bits and pieces (actionItems, views) know
## about their inter-relations, it is only the tree that links the
## individual views and the common evaluation environment containing all
## stored objects. Nodes in the tree correspond to views, edges to
## actionItems, stored as references in the edgeData slot. The tree itself
## is also stored in the environment, which allows for reference-based
## updating without the necessaty of an assignment method or the like.
## In addition to the tree we store an alias table in the environment.
## Internally, all objects are referenced by their guid, but we allow for
## more human readble aliases (usually the "name" slot) which can ge used
## to identify objects if they are unique. Whenever possible we try to plot
## these readable names and those are also available for completion.
## Note that the environment in the prototype gets created once and all
## objects created via "new" without explicitely defining "env" will
## essentially share a common environment. This is fixed in the constructor
## ---------------------------------------------------------------------------
setClass("workFlow",
         representation=representation(name="character",
         tree="fcTreeReference",
         alias="fcAliasReference",
         env="environment"),
         prototype=prototype(name="default",
         tree=fcNullReference(),
         alias=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor takes a flow data object (flowFrame or flowSet) and
## makes a copy in the evaluation environment. It also sets up the views
## graph and the alias table in the environment
workFlow <- function(data, name="default", env=new.env(parent=emptyenv()))
{
    if(!is(data, "flowFrame") && !is(data, "flowSet"))
        stop("'data' must be a flow data structure (flowFrame or flowSet)",
             call.=FALSE)
    ## some sanity checks up front
    checkClass(name, "character", 1)
    checkClass(env, "environment")
    wf <-  new("workFlow", name=name, env=env)
    ## set up the alias table as an environment in the workFlow
    aliasTable <- new.env(hash=TRUE, parent=emptyenv())
    id <- refName(aliasTable)
    assign("alias", id, aliasTable)
    assign(id, aliasTable, wf@env)
    wf@alias <- new("fcAliasReference", env=wf@env, ID=id)
    ## Assign the data to the workFlow and create a base view
    dataRef <- assign(value=data, envir=wf)
    viewRef <- view(workflow=wf, name="base view", data=dataRef)
    ## Set up the views tree
    tree <- new("graphNEL", nodes=identifier(viewRef), edgemode="directed")
    wf@tree <- assign(value=tree, envir=wf)
    return(wf)
}


## Create a reference by assigning 'value' to the symbol 'x' in the
## evaluation environment in 'workflow'. If 'x' is not specified, we
## create a reasonable unique identifier. These methods return a reference
## to 'value' for further use.
## Note that creation of a NULL reference does not result in any assignment
## to the environment, but still a fcNullReference object is returned.


## Assign any object to a workFlow object, the symbol (as guid) is created
## automatically
setMethod("assign",
          signature=signature(x="missing",
          value="ANY",
          pos="workFlow",
          envir="missing",
          inherits="missing",
          immediate="missing"),
          definition=function(value, pos)
      {
          id <- refName(value)
          if(!is.null(value))
              assign(id, value, envir=pos)
          a <- do.call(refType(value), list(ID=id, env=pos))
          return(a)
      })

## The same behaviour as above, but allow workflow to be the 'envir' argument
setMethod("assign",
          signature=signature(x="missing",
          value="ANY",
          pos="missing",
          envir="workFlow",
          inherits="missing",
          immediate="missing"),
          definition=function(value, envir) assign(value=value, pos=envir))

## Assign to a particular symbol (potentially overwriting existing ones)
setMethod("assign",
          signature=signature(x="character",
          value="ANY",
          pos="workFlow",
          envir="missing",
          inherits="missing",
          immediate="missing"),
          definition=function(x, value, pos)
      {
          rmAlias(x, pos)
          if(!is.null(value)){
              if(x %in% ls(pos))
                  warning("Overwriting object in the environment.", call.=FALSE)
              assign(x, value, envir=pos@env)
          }else{
              rm(list=x, envir=pos@env)
          }
          do.call(refType(value), list(ID=x, env=pos))     
      })

## Assign via existing reference.
setMethod("assign",
          signature=signature(x="fcReference",
          value="ANY",
          pos="workFlow",
          envir="missing",
          inherits="missing",
          immediate="missing"),
          definition=function(x, value, pos)
      {
          rmAlias(identifier(x), pos)
          if(is.null(value)){
              Rm(x, rmRef=FALSE)
          }else{
              assign(identifier(x), value, envir=pos@env)  
          }
          do.call(refType(value), list(ID=identifier(x), env=pos))
      })

## The same behaviour as above, but allow workflow to be the 'envir' argument
setMethod("assign",
          signature=signature(x="ANY",
          value="ANY",
          pos="missing",
          envir="workFlow",
          inherits="missing",
          immediate="missing"),
          definition=function(x, value, envir)
      {
          assign(x=x, value=value, pos=envir)
      })



## ===========================================================================
## actionItem
## ---------------------------------------------------------------------------
## actionItems are either gates, transformations or compensations that
## work on a particular view of the data, linked to it by the parentView
## slot.
## ---------------------------------------------------------------------------
setClass("actionItem", 
         representation=representation("VIRTUAL",
         ID="character",
         name="character",
         parentView="fcViewReference",
         alias="fcAliasReference",
         env="environment"),
         prototype=prototype(ID=paste("actionRef", guid(), sep="_"),
         name="",
         parentView=fcNullReference(),
         alias=fcNullReference(),
         env=new.env(parent=emptyenv())))



## ===========================================================================
## gateActionItem
## ---------------------------------------------------------------------------
## A subclass od actionItem. This contains the definiton of the filter/gate
## operation
## ---------------------------------------------------------------------------
setClass("gateActionItem",
         contains="actionItem",
         representation=representation(gate="fcFilterReference",
         filterResult="fcFilterResultReference"),
         prototype=prototype(ID=paste("gateActionRef", guid(), sep="_"),
         gate=fcNullReference(),
         filterResult=fcNullReference()))

## The constructor creates the gateActionItem object and directly assigns
## it to the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
gateActionItem <- function(ID=paste("gateActionRef", guid(), sep="_"),
                           name=paste("action", identifier(get(gate)), sep="_"),
                           parentView, gate, filterResult, workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(gate, "fcFilterReference")
    checkClass(parentView, "fcViewReference")
    if(missing(filterResult))
        filterResult <-fcNullReference()
    action <- new("gateActionItem", ID=ID, name=name, gate=gate,
                  parentView=parentView, env=workflow@env,
                  filterResult=filterResult, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}



## ===========================================================================
## transformActionItem
## ---------------------------------------------------------------------------
## transformation actionItem. This contains the definiton of the
## transformation, which is not defined yet. For now, we use the
## transformMapList class for this purpose
## ---------------------------------------------------------------------------
setClass("transformActionItem",
         contains="actionItem",
         representation=representation(transform="fcTransformReference"))

## The constructor creates the transformActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
transformActionItem <- function(ID=paste("transActionRef", guid(), sep="_"),
                                name="no scheme yet", parentView, transform,
                                workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(transform, "fcTransformReference")
    checkClass(parentView, "fcViewReference")
    action <- new("transformActionItem", ID=ID, name=name,
                  transform=transform, parentView=parentView,
                  env=workflow@env, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}


## ===========================================================================
## compensateActionItem
## ---------------------------------------------------------------------------
## compensation actionItem. This contains the definiton of the
## compensation
## ---------------------------------------------------------------------------
setClass("compensateActionItem",
         contains="actionItem",
         representation=representation(compensate="fcCompensateReference"))

## The constructor creates the compensateActionItem object and directly
## assigns it to the evaluation ennvironment in 'workflow'. The return
## value is a reference to that object.
compensateActionItem <- function(ID=paste("compActionRef", guid(), sep="_"),
                                 name=paste("action", identifier(get(compensate)),
                                 sep="_"),
                                 parentView, compensate,
                                 workflow)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(compensate, "fcCompensateReference")
    checkClass(parentView, "fcViewReference")
    action <- new("compensateActionItem", ID=ID, name=name,
                  compensate=compensate, parentView=parentView,
                  env=workflow@env, alias=workflow@alias)
    return(assign(x=ID, value=action, envir=workflow))
}

## ===========================================================================
## view
## ---------------------------------------------------------------------------
## The concept of views in this context does not adhere strictly to the
## definition. Whenever possible, we try to provide a "real" view on the
## original (or at least on the parent) data, however operations like
## compensation or transformation will alter the data values, and unless
## we want to add new columns to the original data (or make a copy while
## retaining the original set), there is no simple way around that. Views
## mainly exist in order to provide a uniform organisational structure
## for the data after an actionItem has been applied. The first node in
## the workflow tree is always a view with a NULL action reference, and
## a non-NULL data reference.
## ---------------------------------------------------------------------------
setClass("view", 
         representation=representation(ID="character",
         name="character",
         action="fcActionReference",
         env="environment",
         alias="fcAliasReference",
         data="fcDataReference"),
         prototype=prototype(ID=paste("view", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         alias=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the view object and directly assigns it to
## the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
view <- function(workflow, ID=paste("viewRef", guid(), sep="_"),
                 name="default", data, action)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    if(missing(data))
        data <- fcNullReference()
    if(missing(action))
        action <- fcNullReference()
    bv <-  new("view", ID=ID, name=name, env=workflow@env,
               action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}



## ===========================================================================
## gateView
## ---------------------------------------------------------------------------
## Gate views store indices of the gating result for further subsetting as a
## list, where each list item contains indices for a single flowFrame. No
## subset will be produced unless another actionItem is applied to the view.
## We need some form of special treatment for gates that produce multiple
## populations. A gateView will always capture the result for only a single
## sub-population, however, the whole filterResult is necessary for plotting
## and a reference to that and the index in the filterResult (i.e., the
## subpopulation) will be stored along with the view . 
## ---------------------------------------------------------------------------
setClass("gateView",
         contains="view",
         representation=representation(indices="list",
         filterResult="fcFilterResultReference",
         frEntry="character"),
         prototype=prototype(ID=paste("gateViewRef", guid(), sep="_"),
         name="",
         filterResult=fcNullReference(),
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the gateView object and directly assigns it to
## the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
gateView <- function(workflow, ID=paste("gateViewRef", guid(), sep="_"),
                     name="default", action, data, indices, 
                     filterResult, frEntry)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(frEntry, "character", 1)
    checkClass(indices, "list")
    checkClass(action, "fcActionReference")
    checkClass(filterResult, "fcFilterResultReference")
    if(missing(data))
        data <- fcNullReference()
    bv <- new("gateView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, indices=indices,
              filterResult=filterResult, frEntry=frEntry,
              alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## check if the parentView of an actionItem object is the result of a
## filtering operation and, if that is the case, that the filter has been
## applied for subsetting (i.e., a new data set has been created)
applyParentFilter <- function(parent, workflow)
{
    pview <- get(parent)
    dataRef <- Data(pview)
    if(is.null(dataRef)){
        parentData <- Data(parent(pview))
        if(is(parentData, "flowFrame"))
            newData <- parentData[pview@indices[[1]],]
        else{
            newData <- parentData[1:length(parentData)]
            for(i in seq_along(pview@indices))
                newData[[i]] <- newData[[i]][pview@indices[[i]],]
        }
        newDataRef <- assign(value=newData, envir=workflow)
        pview@data <- newDataRef
        assign(parent, value=pview, envir=workflow)
    }
}

## constructor directly from a filter object. This creates a gateActionItem
## in the workFlow and from that a gateView which is also directly stored in
## the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="concreteFilter"),
          definition=function(wf, action, parent=NULL)
      {
          if(is(action, "filterResult"))
              stop("Don't know how to handle object of class '",
                   class(action), "'", call.=FALSE)
          else if(is(action, "filter")){
              ## assign the filter to the evaluation environment and create
              ## a reference to it
              gateRef <- assign(value=action, envir=wf)
              ## get the parentView. If not explicitely specified, use the
              ## root node
              pid <- if(is.null(parent)) views(wf)[[1]] else parent
              pid <- getAlias(pid, wf)
              pview <- fcViewReference(ID=pid, env=wf)
              ## create and assign a new gateActionItem
              actionRef <- gateActionItem(parentView=pview, gate=gateRef,
                                          workflow=wf)
              ## check if the previous filter has been applied for subsetting
              applyParentFilter(pview, wf)
              ## now evaluate the filter and assign the result
              fres <- filter(Data(get(pview)), action)
              if(!validFilterResultList(fres))
                  stop("Don't know how to proceed.", call.=FALSE)
              fresRef <-  assign(value=fres, envir=wf)
              gAction <- get(actionRef)
              gAction@filterResult <- fresRef
              assign(actionRef, gAction, wf) 
              ## we need to distinguish between logicalFilterResults and
              ## multipleFilterResults
              nodes <- NULL
              if(!is(fres, "filterResultList")){
                  len <-
                      if(is(fres, "logicalFilterResult")) 2
                      else length(fres)
                  for(i in seq_len(len)){
                      vid <- gateView(name=names(fres)[i], workflow=wf,
                                      action=actionRef, filterResult=fresRef,
                                      indices=list(fres[[i]]@subSet),
                                      frEntry=names(fres)[i])
                      nodes <- c(nodes, identifier(vid))
                  }
              }else{
                  len <-
                      if(is(fres[[1]], "logicalFilterResult")) 2
                      else length(fres[[1]])
                  for(i in seq_len(len)){
                      vid <- gateView(name=names(fres[[1]])[i], workflow=wf,
                                      action=actionRef, filterResult=fresRef,
                                      indices=lapply(fres, function(y) y[[i]]@subSet),
                                      frEntry=names(fres[[1]])[i])
                      nodes <- c(nodes, identifier(vid))
                  }
              }
              ## update the filter and filterResult IDs
              identifier(action) <- paste("filter", identifier(action),
                                          sep="_")
              identifier(fres) <- paste("fres", identifier(fres),
                                        sep="_")
              assign(gateRef, action, wf)
              assign(fresRef, fres, wf)
              ## add new nodes and edges to the workflow tree
              tree <- get(wf@tree)
              tree <- addNode(nodes, tree)
              tree <- addEdge(pview@ID, nodes, tree)
              edgeDataDefaults(tree, "actionItem") <- fcNullReference()
              edgeData(tree, pview@ID, nodes, "actionItem") <-
                  actionRef
              assign(x=wf@tree, value=tree, envir=wf)
              return(wf)
          } else stop("Don't know how to handle object of class '",
                      class(action), "'", call.=FALSE)       
      })





## ===========================================================================
## transformView
## ---------------------------------------------------------------------------
## A transfomation makes a copy of the data independent of whether it
## is a leaf node or not.
## Do we want to allow to introduce new data columns. My guess is that we
## have to, but for now let's assume we don't.
## ---------------------------------------------------------------------------
setClass("transformView",
         contains="view",
         prototype=prototype(ID=paste("transViewRef", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the transformView object and directly assigns it to
## the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
transformView <- function(workflow, ID=paste("transViewRef", guid(), sep="_"),
                          name="default", action, data)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(action, "fcActionReference")
    checkClass(data, "fcDataReference")
    bv <- new("transformView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## constructor directly from a transformation object. This creates a
## transformActionItem in the workFlow and from that a transformView
## which is also directly stored in the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="transformList"),
          definition=function(wf, action, parent=NULL)
      {
          ## assign the transformation to the evaluation environment and
          ## create a reference to it
          transRef <- assign(value=action, envir=wf)
          ## get the parentView. If not explicitely specified, use the
          ## root node
          pid <- if(is.null(parent)) views(wf)[[1]] else parent
          pid <- getAlias(pid, wf)
          pview <- fcViewReference(ID=pid, env=wf)
          tree <- get(wf@tree)
          if(length(unlist(adj(tree, pid))))
              warning("The selected parent view is not a leaf node.\n",
                      "Don't know how to update yet.", call.=FALSE)
          ## create and assign a new transformActionItem
          actionRef <- transformActionItem(parentView=pview,
                                           transform=transRef,
                                           workflow=wf)
          ## check if the previous filter has been applied for subsetting
          applyParentFilter(pview, wf)
          ## now transform the data and assign the result
          tData <- action %on% Data(get(pview))
          dataRef <- assign(value=tData, envir=wf)
          vid <- transformView(name="no scheme yet", workflow=wf,
                               action=actionRef, data=dataRef)
          
          ## add new nodes and edges to the workflow tree
          nid <- identifier(vid)
          tree <- addNode(nid, tree)
          tree <- addEdge(pid, identifier(vid), tree)
          edgeDataDefaults(tree, "actionItem") <- fcNullReference()
          edgeData(tree, pid , nid, "actionItem") <- actionRef
          assign(x=wf@tree, value=tree, envir=wf)
          return(wf)   
      })



## ===========================================================================
## compensateView
## ---------------------------------------------------------------------------
## A compensation makes a copy of the data independent of whether it
## is a leaf node or not.
## ---------------------------------------------------------------------------
setClass("compensateView",
         contains="view",
         prototype=prototype(ID=paste("compViewRef", guid(), sep="_"),
         name="",
         action=fcNullReference(),
         data=fcNullReference(),
         env=new.env(parent=emptyenv())))

## The constructor creates the compensateView object and directly assigns it
## to the evaluation ennvironment in 'workflow'. The return value is a
## reference to that object.
compensateView <- function(workflow, ID=paste("compViewRef", guid(), sep="_"),
                           name="default", action, data)
{
    checkClass(workflow, "workFlow")
    checkClass(ID, "character", 1)
    checkClass(name, "character", 1)
    checkClass(action, "fcActionReference")
    checkClass(data, "fcDataReference")
    bv <- new("compensateView", ID=ID, name=name, env=workflow@env,
              action=action, data=data, alias=workflow@alias)
    ref <- assign(identifier(bv), bv, workflow)
    return(ref)
}

## constructor directly from a compensation object. This creates a
## compensateActionItem in the workFlow and from that a compensateView
## which is also directly stored in the workFlow object.
setMethod("add",
          signature=signature(wf="workFlow", action="compensation"),
          definition=function(wf, action, parent=NULL)
      {
          ## assign the compensation to the evaluation environment and
          ## create a reference to it
          compRef <- assign(value=action, envir=wf)
          ## get the parentView. If not explicitely specified, use the
          ## root node
          pid <- if(is.null(parent)) views(wf)[[1]] else parent
          pid <- getAlias(pid, wf)
          if(pid != getAlias(views(wf), wf)[1])
              warning("The selected parent view is not a root node.\n",
                      "Are you sure this is correct?", call.=FALSE)
          pview <- fcViewReference(ID=pid, env=wf)
          
          ## create and assign a new ActionItem
          actionRef <- compensateActionItem(parentView=pview,
                                            compensate=compRef,
                                            workflow=wf)      
          ## check if the previous filter has been applied for subsetting
          applyParentFilter(pview, wf)
          ## now transform the data and assign the result
          tData <- compensate(Data(get(pview)), action)
          dataRef <- assign(value=tData, envir=wf)
          vid <- compensateView(name=identifier(action), workflow=wf,
                                action=actionRef, data=dataRef)
          ## update the identifier of the compensation object
          identifier(action) <- paste("comp", identifier(action),
                                      sep="_")
          assign(compRef, value=action, envir=wf) 
          ## add new nodes and edges to the workflow tree
          nid <- identifier(vid)
          tree <- get(wf@tree)
          tree <- addNode(nid, tree)
          tree <- addEdge(pid, identifier(vid), tree)
          edgeDataDefaults(tree, "actionItem") <- fcNullReference()
          edgeData(tree, pid , nid, "actionItem") <- actionRef
          assign(x=wf@tree, value=tree, envir=wf)
          return(wf)   
      })


## ===========================================================================
## Unity transformation
## ---------------------------------------------------------------------------
## Transforms parameters names provided as characters into unity transform 
## objects which can be evaluated to retrive the corresponding columns from the
## data frame
## ---------------------------------------------------------------------------
setClass("unitytransform",
	 contains="transform",
	 representation=representation(parameters="character"))

unitytransform <- function(parameters,
                           transformationId="defaultUnityTransform")
{
    checkClass(transformationId, "character", 1)
    if(missing(parameters))
        parameters <- character()
    new("unitytransform", parameters=parameters,
        transformationId=transformationId)
}



## ===========================================================================
## Polynomial transformation of degree 1 
## ---------------------------------------------------------------------------
## Allows for scaling ,linear combination and translation within a single 
## transformation
## ---------------------------------------------------------------------------
setClass("dg1polynomial", 		
         contains="transform",
         representation=representation(parameters="parameters",
                                       a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=new("parameters"),
                             a=1,
                             b=1))

dg1polynomial <- function(parameters, a=1, b=1,
                          transformationId="defaultDg1polynomialTransform")
{
    checkClass(a, "numeric", length(parameters))
    checkClass(b, "numeric", 1)
    checkClass(transformationId, "character", 1)
    new("dg1polynomial", parameters=parameters, a=a, b=b,
        transformationId=transformationId)
}



## ===========================================================================
## Ratio transformation
## ---------------------------------------------------------------------------
## Ratio of two arguments defined in the transformation
## ---------------------------------------------------------------------------
setClass("ratio",
         contains="transform",
         representation(numerator="transformation",
                        denominator="transformation"),
	 prototype=prototype(numerator=unitytransform(),
                             denominator=unitytransform()))

ratio <- function(numerator=unitytransform(),
                  denominator=unitytransform(),
                  transformationId="defaultRatioTransform")
{
    if(!is(numerator, "transform")){
        checkClass(numerator, "character", 1)
        numerator <- unitytransform(numerator)
    }
    if(!is(denominator, "transform")){
        checkClass(denominator, "character", 1)
        denominator=unitytransform(denominator)
    }  
    new("ratio", numerator=numerator, denominator=denominator,
        transformationId=transformationId)
}



## ===========================================================================
## Quadratic transformation
## ---------------------------------------------------------------------------
setClass("quadratic", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             a=1),
         validity=function(object) 
     {
         msg<-NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Quadratic transform is defined for one parameter")
         if(length(object@a)!=1)
             msg <- c(msg, "Only one coefficient is defined for quadratic transform")
         if(object@a==0)
             msg <- c(msg, "'a' should be non zero")
         msg
     })

quadratic <- function(parameters="NULL", a=1,
                      transformationId="defaultQuadraticTransform")
    new("quadratic",parameters=parameters,a=a,
        transformationId=transformationId)

          

## ===========================================================================
## Squareroot transformation
## ---------------------------------------------------------------------------
setClass("squareroot", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Square root transform is defined for one parameter")
         if(length(object@a)!=1)
             msg <- c(msg, "Only one coefficient is defined for quadratic transform")
         if(object@a==0)
             msg <- c(msg, "Coefficient should be non zero")
         msg
     })

squareroot <- function(parameters, a=1,
                       transformationId="defaultSquarerootTransform")
    new("squareroot", parameters=parameters, a=a,
        transformationId=transformationId)



## ===========================================================================
##  Loarithmical Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
setClass("logarithm",
         contains="singleParameterTransform",
         representation=representation(a="numeric", b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Logarithm transform is defined for one parameter")
         if(object@a==0)
             msg <- c(msg, "'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non zero number")
         msg
     })

logarithm <- function(parameters, a=1, b=1,
                      transformationId="defaultLogarithmTransform")
    new("logarithm", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Exponential Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------

setClass("exponential", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             a=1,
                             b=1),
         validity=function(object) 
     {
         msg <-NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Exponential transform is defined for one parameter")
         if(object@a==0)
             msg<-c(msg,"'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg,"'b' should be a non zero number")
         msg  
     })

exponential <- function(parameters, a=1, b=1,
                        transformationId="defaultExponentialTransformation")
    new("exponential", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Inverse hyperbolic sin Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
setClass("asinht", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Inverse hypberbolic transform is defined for one parameter")
         if(object@a==0)
             msg <- c(msg, "'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non zero number")
         msg
     })

asinht <- function(parameters="NULL", a=1, b=1,
                   transformationId="defaultAsinhTransform")
    new("asinht", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Inverse hyperbolic sin Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## ---------------------------------------------------------------------------
setClass("sinht", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             a=1,
                             b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Hypberbolic transform is defined for one parameter")
         if(object@a==0)
             msg <- c(msg, "'a' should be a non zero number")
         if(object@b==0)
             msg <- c(msg, "'b' should be a non zero number")
         msg
     })

sinht <- function(parameters, a=1, b=1,
                  transformationId="defaultSinhtTransform")
    new("sinht", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  Hyperlog Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## --------------------------------------------------------------------------- 
setClass("hyperlog", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Hyperlog transform is defined for one parameter")
         if(object@a<=0)
             msg <- c(msg, "'a' should be greater than zero")
         if(object@b<=0)
             msg <- c(msg, "'b' should be greater than zero")
         msg
     })

hyperlog <- function(parameters="NULL", a=1, b=1,
                     transformationId="defaultHyperlogTransform")
    new("hyperlog", parameters=parameters, a=a, b=b,
        transformationId=transformationId)



## ===========================================================================
##  EH Transformation 
## ---------------------------------------------------------------------------
## inputs a,b of type numeric and parameter of type transformation or character
## --------------------------------------------------------------------------- 
setClass("EHtrans", 		
         contains="singleParameterTransform",
         representation=representation(a="numeric",
                                       b="numeric"),
         prototype=prototype(parameters=unitytransform(), a=1, b=1),
         validity=function(object) 
     {
         msg <-NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "EH transform is defined for one parameter")
         if(object@a<=0)
             msg<-c(msg, "'a' should be greater than zero")
         if(object@b<=0)
             msg<-c( msg, "'b' should be greater than zero")
         msg
     })

EHtrans <- function(parameters, a=1, b=1,
                    transformationId="defaultEHtransTransform")
    new("EHtrans", parameters=parameters, a=a, b=b,
        transformationId=transformationId)
      
          

## ===========================================================================
##  Splitscale Transformation 
## ---------------------------------------------------------------------------
setClass("splitscale", 		
         contains="singleParameterTransform",
         representation=representation(r="numeric",
                                       maxValue="numeric",
                                       transitionChannel="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             r=1,
                             maxValue=1,
                             transitionChannel=4),
         validity=function(object) 
     {
         msg <-NULL
         if(length(object@parameters)!=1)
             msg<-c(msg, "Split scale transform is defined for one parameter")
         if(object@r<=0)
             msg <- c(msg, "'r' should be a greater than zero")
         if(object@maxValue<=0)
             msg <- c(msg, "maxValue should be a greater than zero")
         if(object@transitionChannel<0)
             msg <- c(msg, "transitionChannel should be a non negative")
         msg
     })

splitscale <- function(parameters="NULL", r=1, maxValue=1, transitionChannel=4,
                       transformationId="defaultSplitscaleTransform")
    new("splitscale",
        parameters=parameters, r=r, maxValue=maxValue,
        transitionChannel=transitionChannel,
        transformationId=transformationId)



## ===========================================================================
##  Inverse Splitscale Transformation 
## ---------------------------------------------------------------------------
setClass("invsplitscale", 		
         contains="singleParameterTransform",
         representation=representation(r="numeric",
                                       maxValue="numeric",
                                       transitionChannel="numeric"),
         prototype=prototype(parameters=unitytransform(),
                             r=1,
                             maxValue=1,
                             transitionChannel=4),
         validity=function(object) 
     {
         msg <- NULL
         if(length(object@parameters)!=1)
             msg <- c(msg, "Split scale transform is defined for one parameter")
         if(object@r<=0)
             msg <- c(msg, "'r' should be a greater than zero")
         if(object@maxValue<=0)
             msg <- c(msg, "'maxValue' should be a greater than zero")
         if(object@transitionChannel<0)
             msg <- c(msg, "'transitionChannel' should be a non negative")
         msg
     })
       
invsplitscale <- function(parameters, r=1, maxValue=1,
                          transitionChannel=4,
                          transformationId="defaultInvsplitscaleTransforms")
    new("invsplitscale",
        parameters=parameters, r=r, maxValue=maxValue,
        transitionChannel=transitionChannel,
        transformationId=transformationId)
      


## ===========================================================================
## Transformation reference
## ---------------------------------------------------------------------------
## Reference to a transformation defined previously
## ---------------------------------------------------------------------------
setClass("transformReference",
         contains="transform",
         representation(searchEnv="environment"))

transformReference <- function(referenceId="defaultTransformReference",
                               searchEnv)
    new("transformReference",
        transformationId=referenceId, searchEnv=searchEnv)
    
