##===========================================================================
## Eval methods for evaluating transformation objects 
## ---------------------------------------------------------------------------


## Helper function to fully resolve all transformationReferences
resolve <- function(x, df)
{
    if(!is(x, "transformReference"))
        eval(x)(df) else resolveTransformReference(x, df)  
}



## ===========================================================================
## Unity transformation
## ---------------------------------------------------------------------------
setMethod("eval",
          signature=signature(expr="unitytransform",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr,envir,enclos)
      {        
          function(df)
              return(df[, expr@parameters, drop=FALSE])
          
      })


## ===========================================================================
## Polynomial transformation of degree 1
## ---------------------------------------------------------------------------

setMethod("eval", 
	  signature=signature(expr="dg1polynomial",
                              envir="missing",
                              enclos="missing"),
	  definition=function(expr, envir, enclos)
      {    
          function(df)
            {   
                par <- expr@parameters
                temp <- sapply(par,function(i){
                    if(is(i, "function")) return(i)
                    resolve(i, df)
                })
                pp <- matrix(expr@a, ncol=length(expr@parameters))
                rowSums(t(apply(temp, 1, function(x) pp*x)))+expr@b       
            }
      })


## ===========================================================================
## Ratio transformation of two arguments
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="ratio",
                              envir="missing",
                              enclos="missing"),
          function(expr,envir,enclos)
      { 
          function(df)
          {  
              num <- resolve(expr@numerator, df)
              den <- resolve(expr@denominator, df)
              num/den
          }
      })


## ===========================================================================
## Quadratic transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
	  signature=signature(expr="quadratic",
                              envir="missing",
                              enclos="missing"),
	  definition=function(expr,envir,enclos)
      {    
          function(df)
          {      
              parameter <- resolve(expr@parameters, df)
              expr@a*(parameter^2)
          }
      })


## ===========================================================================
## Squareroot transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
	  signature=signature(expr="squareroot",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr,envir,enclos)
      {
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              sqrt(abs(parameter/(expr@a)))
          }
      })


## ===========================================================================
## Log transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
	  signature=signature(expr="logarithm",
                              envir="missing",
                              enclos="missing"),
	  definition=function(expr, envir, enclos)
      {    
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              temp <- expr@a*parameter
              result <- vector(mode="numeric", length=length(temp))
              result[temp>0] <- log(temp[temp>0])*expr@
              return(result)
          }
      })  


## ===========================================================================
## Exponential transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="exponential",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr, envir, enclos)
      {
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              exp(parameter/expr@b)/expr@a
          }
      })


## ===========================================================================
## Inverse hyperbolic sin transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="asinht",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr, envir, enclos)
      {    
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              expr@b*asinh(expr@a*parameter)
          }
      })


## ===========================================================================
## Hyperbolic sin transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="sinht",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr, envir, enclos)
      {
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              sinh(parameter/(expr@a))/expr@b
          }
      }) 


## ===========================================================================
## Hyperlog transformation 
## ---------------------------------------------------------------------------


setMethod("eval", 
	  signature=signature(expr="hyperlog",
                              envir="missing",
                              enclos="missing"),
	  definition=function(expr, envir, enclos)
	  {    
            function(df)
              {
                  parameter <- resolve(expr@parameters, df)
                  solveEH(expr@parameters , expr@a ,expr@b, df)
               }
	  })


## ===========================================================================
## Splitscale transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="splitscale",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr, envir, enclos)
      {    
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              args <- expr@parameters
              transitionChannel <- expr@transitionChannel
              r <- expr@r
              maxValue <- expr@maxValue
              b <- transitionChannel/2
              d <- 2*r*log10(2.71828)/transitionChannel
              log10t <- -2*log10(2.71828)*r/transitionChannel +log10(maxValue)
              t <- 10^(log10t)
              a <- transitionChannel/(2*t)
              log10ct <- (a*t+b)*d/r
              c <- (10^log10ct)/t
              idx <- which(args <= t)
              idx2 <- which(args > t)
              if(length(idx2)>0)
                  args[idx2] <- log10(c*args[idx2])*r/d
              if(length(idx)>0)
                  args[idx] <- a*args[idx]+b
              args
          }
      })


## ===========================================================================
## Inverse Splitscale transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="invsplitscale",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr, envir, enclos)
      {    
          function(df)
          {
              parameter <- resolve(expr@parameters, df)
              args <- expr@parameters
              transitionChannel <- expr@transitionChannel
              r <- expr@r
              maxValue <- expr@maxValue
              b <- transitionChannel/2
              d <- 2*r*log10(2.71828)/transitionChannel
              log10t <- -2*log10(2.71828)*r/transitionChannel +log10(maxValue)
                  t <- 10^(log10t)
                  a <- transitionChannel/(2*t)
                  log10ct <- (a*t+b)*d/r
                  c <- (10^log10ct)/t            
                  thresh <- t*a+b
                  idx <- which(args<=thresh)
                  idx2 <- which(args>thresh)
                  if(length(idx2)>0)
                      args[idx2] <- (10^(args[idx2]*(d/r)))/c
                  if(length(idx>0))
                      args[idx] <- (args[idx]-b)/a
                  args
              }
        })


## ===========================================================================
## Transformation reference
## ---------------------------------------------------------------------------
setMethod("eval",
          signature=signature(expr="transformReference",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr,envir,enclos)
              expr@searchEnv[[slot(expr,"transformationId")]]
          )	
	
####----------------------------------------------------------------------------
#### 
#### Reference to a predefined gate 
#### 
####----------------------------------------------------------------------------

setMethod("eval",
          signature=signature(expr="filterReference",
                              envir="missing",
                              enclos="missing"),
          definition=function(expr,envir,enclos)
          expr@env[[expr@name]]  #retrieved using name instead of the filterId
          )	
