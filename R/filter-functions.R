## ==========================================================================
## function to print filter summary
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print.filterSummary <- function(x,...) {
    
    if(length(x$name) == 1) {
		with(x,cat(sprintf("%s: %d of %d (%.2f%%)\n",name,true,n,100*p)))
	} else {
	for(i in seq(along=x$name))
			with(x,cat(sprintf("%s: %d of %d (%.2f%%)\n",name[i],true[i],n[i],100*p[i])))
	}
}


## ==========================================================================
## fitNorm2
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitNorm2 <- function(x, y=NA, scalefac=1, method="covMcd", noise,
                     gateName="fitNorm") {

  if(!require(rrcov))
    stop("Required package rrcov could not be found.")
  if(is(x, "cytoFrame"))
    x <- exprs(x)[,1:2]
    
  if (!(is.matrix(x) && ncol(x)==2)){
    if (!length(x)==length(y) || !is.numeric(x) || !is.numeric(y))
      stop("'x' and 'y' must be numeric vectors of equal length")
    x <- cbind(x, y)
  }
  xorig <- x
  if (!missing(noise)){
    if(is.logical(noise))
      noise <- which(noise)
    if (!is.numeric(noise) || length(noise) > nrow(x) || length(noise)==0)
      stop("'noise' should be an index or logical vector not longer than x") 
    x <- x[-noise, ,drop=FALSE]
  }
  if (nrow(x)<50)
    stop("Not enough data points for reliable analysis")
  
  if (!is.numeric(scalefac))
    stop("'scalefac' must be numeric")
  
  cov <- switch(method,
    covMcd = {
      nmax <- 50000
      if (nrow(x)>nmax)
        covMcd(x[sample(nrow(x), nmax),])
      else
        covMcd(x)
    },
    cov.rob = {
      cov.rob(x)
    },
    stop("'method' must be one of 'covMcd' or 'cov.rob'")
   ) ## end of switch
  
  mu   <- cov$center
  S    <- cov$cov
  Sinv <- solve(S)
  w    <- rbind(xorig[,1], xorig[,2])-mu
  z    <- Sinv %*% w
  p    <- exp(-0.5 * (z[1,]*w[1,] +  z[2,]*w[2,]))
  sel  <- p > exp(-0.5 * scalefac^2)

  gfun <- function(x=x, cov, scalefac){
    mu   <- cov$center
    S    <- cov$cov
    Sinv <- solve(S)
    w    <- rbind(x[,1], x[,2])-mu
    z    <- Sinv %*% w
    p    <- exp(-0.5 * (z[1,]*w[1,] +  z[2,]*w[2,]))
    return(p > exp(-0.5 * scalefac^2))
  }
  cn <- colnames(x)
  if(is.null(cn))
    colnames(xorig) <- c("x", "y")
     
  gate <- new("gate", name=gateName,
              gateFun=function(x) gfun(x=x, cov=cov, scalefac=scalefac),
              colnames=colnames(xorig),
              logic="&", type="fitNorm") 
  return(invisible(list(mu=mu, S=S, p=p, sel=sel, scalefac=scalefac,
                        data=xorig, gate=gate)))
  
}









