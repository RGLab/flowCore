## ==========================================================================
## Linear transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
linearTransform <- function(transformationId,a=1,b=0){
    if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
    function(x){    
        x <- a*x+b
    }
}

## ==========================================================================
## Quadratic transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
quadraticTransform <- function(transformationId,a,b,c){
  if(!is.double(a)) 
      stop("a must be numeric")
    if(!is.double(b))
       stop("b must be numeric")
  if(!is.double(c))
       stop("c must be numeric")
    function(x){
        x <- a*x^2 + b*x + c
    }
}

## ==========================================================================
## Natural logarithm transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lnTransform <- function(transformationId,r,d){
    if(!is.double(r) || r <= 0)
       stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
       stop("d must be numeric")
    function(x){
     x<-log(x)*(r/d)
 }
}

## ==========================================================================
## Logarithm transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
logTransform <- function(transformationId,logbase=10,r,d){
     if(!is.double(r) || r <= 0)
       stop("r must be numeric and positive")
    if(!is.double(d) || d <=0)
       stop("d must be numeric")
      if(!is.double(r) || r <=0)
       stop("r must be numeric and positive")
    if(!is.double(logbase) || logbase <= 1)
       stop("logabse must be a pnumeric greater than 1")
    function(x){
        x <- log(x,logbase)*(r/d)
    }
}

## ==========================================================================
## General biexponential transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
biexponentialTransform<- function(transformationId,a=.5,b=1,c=.5,d=1,f=0,w=0,
           tol=.Machine$double.eps^0.25,maxit=as.integer(5000)){
    function(x){
        x <- .Call(biexponential_transform,x,a,b,c,d,f,w,tol,maxit)
    }
}

logicleTransform <- function(w=0,r=262144,d=5,...) {
  if(w>d) stop("Negative range decades must be smaller than total number of decades")
  w = w*log(10)
  d = d*log(10)
  p = if(w==0) 1 else uniroot(function(p) -w+2*p*log(p)/(p+1),c(.Machine$double.eps,2*(w+d)))$root
  ##new("biexponentialTransformation",a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...)
 biexponentialTransform(a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...)
}

## ===========================================================================
## Truncation Transformation
## ---------------------------------------------------------------------------
truncateTransform <- function(transformationId,a){
    function(x){
        x[x<a] <- a
        x
    }
}


## ===========================================================================
## Scale Transformation
## ---------------------------------------------------------------------------
scaleTransform <- function(transformationId,a,b){
    function(x){
     	x=(x-a)/(b-a)
    }
}

## ===========================================================================
## Hyperbolic Arcsin Transformation
## ---------------------------------------------------------------------------
arcsinhTransform <- function(transformationId,a,b,c=1) {
	function(x) asinh(a+b*x)+log(c)
}

