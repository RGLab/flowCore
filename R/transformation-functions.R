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
    function(x){
        x <- a*x^2 + b*x + c
    }
}

## ==========================================================================
## Natural logarithm transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lnTransform <- function(transformationId,r,d){
    function(x){
     x<-log(x)*(r/d)
 }
}

## ==========================================================================
## Logarithm transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
logTransform <- function(transformationId,logbase,r,d){
    function(x){
        x <- log(x,logbase)*(r/d)
    }
}

## ==========================================================================
## General biexponential transformation function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
biexponentialTransform<- function(transformationId,a,b,c,d,f,w,tol,maxit){
    function(x){
        x <- .Call(biexponential_transform,x,a,b,c,d,f,w,tol,maxit)
    }
}

logicleTransform <- function(w=0,r=262144,d=5,...) {
  if(w>d) stop("Negative range decades must be smaller than total number of decades")
  w = w*log(10)
  d = d*log(10)
  p = if(w==0) 1 else uniroot(function(p) -w+2*p*log(p)/(p+1),c(.Machine$double.eps,2*(w+d)))$root
  new("biexponentialTransformation",a=r*exp(-(d-w)),b=1,c=r*exp(-(d-w))*p^2,d=1/p,f=p^2-1,w=w,...)
}
