## ==========================================================================
## Linear transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
linearTransformation <- function(transformation,flowObject){
    parameters=transformation@parameters
    a=transformation@a
    b=transformation@b
    ## test for validity of objects
    if(class(flowObject)!="flowFrame"){
        stop("flowObject must of class flowFrame.")
    }
    if (is.null(flowObject@exprs)) {
        stop("There is no data to createFilter.")
    }
    data <- flowObject@exprs
    flowObject@exprs[,parameters] <- a * data[,parameters] + b
    flowObject@description[["transformation"]] <- paste("A linear transformation was applied on the",
                                                     parameters, "channel",sep="",collapse=" ")
    
    return(flowObject)
}

## ==========================================================================
## Quadratic transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
quadraticTransformation <- function(transformation,flowObject){
    parameters=transformation@parameters
    a=transformation@a
    b=transformation@b
    c=transformation@c
    ## test for validity of objects
    
    if(class(flowObject)!="flowFrame"){
        stop("flowObject must of class flowFrame.")
    }
    if (is.null(flowObject@exprs)) {
        stop("There is no data to createFilter.")
    }
    data <- flowObject@exprs
    flowObject@exprs[,parameters] <- a * data[,parameters]^2 + b * data[,parameters] + c
    flowObject@description[["transformation"]] <- paste("A quadratic transformation was applied on the",
                                                     parameters, "channel",sep="",collapse=" ")
    
    return(flowObject)
}

## ==========================================================================
## Natural logarithm transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lnTransformation <- function(transformation,flowObject){
    parameters=transformation@parameters
    r=transformation@r
    d=transformation@d
   
    ## test for validity of objects
    if(class(flowObject)!="flowFrame"){
        stop("flowObject must of class flowFrame.")
    }
    if (is.null(flowObject@exprs)) {
        stop("There is no data to createFilter.")
    }

    data <- flowObject@exprs
    flowObject@exprs[,parameters] <- log(data[,parameters])*(r/d)
    flowObject@description[["transformation"]] <- paste("A natural logarithm transformation was applied on the",
                                                     parameters, "channel",sep="",collapse=" ")
    
    return(flowObject)
}
    
    
## ==========================================================================
## Logarithm transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
logTransformation <- function(transformation,flowObject){
    parameters=transformation@parameters
    r=transformation@r
    d=transformation@d
    logbase=transformation@logbase
   
    ## test for validity of objects
    if(class(flowObject)!="flowFrame"){
        stop("flowObject must of class flowFrame.")
    }
    if (is.null(flowObject@exprs)) {
        stop("There is no data to createFilter.")
    }
    
    data <- flowObject@exprs
    flowObject@exprs[,parameters] <- log(data[,parameters],base=logbase)*(r/d)
    flowObject@description[["transformation"]] <- paste("A log transformation was applied on the",
                                                     c(parameters), "channel",sep="",collapse=" ")
    
    return(flowObject)
}
    

## ==========================================================================
## General biexponential transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
biexponentialTransformation <- function(transformation,flowObject) {
  parameters=transformation@parameters
  a=transformation@a
  b=transformation@b
  c=transformation@c
  d=transformation@d
  f=transformation@f
  w=transformation@w
  tol=transformation@tol
  maxit=transformation@maxit
  
  ## test for validity of objects
  if(class(flowObject)!="flowFrame"){
    stop("flowObject must of class flowFrame.")
  }
  if (is.null(flowObject@exprs)) {
    stop("There is no data to createFilter.")
  }
  flowObject@exprs[,parameters] <- .Call(biexponential_transform,flowObject@exprs[,parameters],a,b,c,d,f,w,tol,maxit)
  flowObject@description[["transformation"]] <- paste("A biexponential transformation was applied on the",
                                                   c(parameters), "channel",sep="",collapse=" ")
  return(flowObject)
}


logicleTransform <- function(w=0,r=262144,d=5,...) {
  if(w>d) stop("Negative range decades must be smaller than total number of decades")
  w = w*log(10)
  d = d*log(10)
  p = if(w==0) 1 else uniroot(function(p) -w+2*p*log(p)/(p+1),c(.Machine$double.eps,2*(w+d)))$root
  new("biexponentialTransformation",a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...)
}
