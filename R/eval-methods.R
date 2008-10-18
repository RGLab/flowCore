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
                coeff <- expr@a
                res <- 0
                for(i in seq_along(expr@parameters))
                  res <- res + temp[,i,drop=FALSE]*coeff[i]
                res <- res+expr@b
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
              result[temp>0] <- log(temp[temp>0])*expr@b
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
              sinh(parameter/(expr@b))/expr@a
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
                  solveEH(parameter , expr@a ,expr@b)
               }
	  })

## ===========================================================================
## EH transformation 
## ---------------------------------------------------------------------------

setMethod("eval", 
	  signature=signature(expr="EHtrans",
                              envir="missing",
                              enclos="missing"),
	  definition=function(expr, envir, enclos)
	  {    
            function(df)
              {  
                  parameter <- resolve(expr@parameters, df)
                  result=0
                  result[parameter>=0]=10^(parameter/expr@a)+ (expr@b*parameter)/expr@a-1
                  result[parameter<0]= -1*(10^(-parameter/expr@a))+ (expr@b*parameter)/expr@a+1
                  return(matrix(result))
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
              idx <- which(parameter <= t)
              idx2 <- which(parameter > t)
              if(length(idx2)>0)
                  parameter[idx2] <- log10(c*parameter[idx2])*r/d
              if(length(idx)>0)
                  parameter[idx] <- a*parameter[idx]+b
              parameter
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
                  idx <- which(parameter<=thresh)
                  idx2 <- which(parameter>thresh)
                  if(length(idx2)>0)
                      parameter[idx2] <- (10^(parameter[idx2]*(d/r)))/c
                  if(length(idx>0))
                      parameter[idx] <- (parameter[idx]-b)/a
                  parameter
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
#### Reference to a predefined gate 
####----------------------------------------------------------------------------

setMethod("eval",
          signature=signature(expr="filterReference",
            envir="missing",
            enclos="missing"),
          definition=function(expr,envir,enclos)
          expr@env[[expr@name]]  #retrieved using name instead of the filterId
          )	


####----------------------------------------------------------------------------
#### Eval transforms for compensations
####----------------------------------------------------------------------------
setMethod("eval",
          signature=signature(expr="compensatedParameter",
            envir="missing",
            enclos="missing"),
          function(expr, envir, enclos)
          {    function(df)
                 { 
                   parameter <- expr@parameters
                   compObj <- expr@searchEnv[[expr@spillRefId]]
                   spillMat <- compObj@spillover
                   cols <- colnames(spillMat)
                   trans <- compObj@parameters
                   for(i in trans){
                     if(is(i, "transformReference")){
                       newCol <- resolveTransformReference(i, df)
                       colnames(newCol) <- i@transformationId
                       df <- cbind(df,newCol)
                     }
                   }
                   t(solve(spillMat)[parameter,]%*%t(df[,cols])) 
                   
                 }
             })

