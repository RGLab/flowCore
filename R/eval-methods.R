##===========================================================================
## Eval methods for evaluating transformation objects 
## ---------------------------------------------------------------------------


## ===========================================================================
## Unity transformation
## ---------------------------------------------------------------------------
setMethod("eval",
          signature=signature(expr="unitytransform",envir="missing",
            enclos="missing"),
          function(expr,envir,enclos)
          {        
              function(df)
                {
                  return(df[,expr@parameters,drop=FALSE])
                }
            }
          )

## ===========================================================================
## Polynomial transformation of degree 1
## ---------------------------------------------------------------------------

setMethod("eval", 
	  signature=signature(expr="dg1polynomial",envir="missing",
            enclos="missing"),
	  function(expr,envir,enclos)
	  {    
            function(df)
	      {   
                result=slot(object,"b")+
                  rowSums(slot(object,"a")*
                          sapply(slot(expr,"parameters"),
                                 function(i)
                                 {   
                                   if(class(i)=="transformReference")
                                     i=eval(i)
                                   eval(i)(df)
                                 }
                                 ) 
                          )                  
              }
	  }
          )

## ===========================================================================
## Ratio transformation of two arguments
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature=signature(expr="ratio",envir="missing",enclos="missing"),
          function(expr,envir,enclos)
          { 
            function(df)
	      {  
                num=slot(expr,"numerator")
                den=slot(expr,"denominator")

                if(class(num)!="transformReference")
                  {
                    num=eval(num)(df)
                  }
                else
                  {
                    num=resolveTransformReference(num,df)
                  }
                if(class(den)!="transformReference")
                  {
                    den=eval(den)(df)
                  }
                else
                  {
                    den=resolveTransformReference(den,df)
                  }    
                result=num/den
              }
	  }
          )
## ===========================================================================
## Quadratic transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
	  signature=signature(expr="quadratic",envir="missing",enclos="missing"),
	  function(expr,envir,enclos)
          {    
            function(df)
              {      
                parameter=slot(expr,"parameters")
                if(class(parameter)!="transformReference")
                  {
                    parameter=eval(parameter)(df)
                  }
                else
                  {
                    parameter=resolveTransformReference(parameter,df)
                  }
                
                result=slot(expr,"a")*(parameter)^2
              }
	  }
          )
## ===========================================================================
## Squareroot transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
	  signature="squareroot",
	  function(expr,envir,enclos)
	  {   function(df)
                {     parameter=slot(expr,"parameters")
                      if(class(parameter)!="transformReference")
                        {
                          parameter=eval(parameter)(df)
                        }
                      else
                        {
                          parameter=resolveTransformReference(parameter,df)
                        } 
                      sqrt(abs(parameter/(slot(expr,"a"))))
                    }
            }
          )

## ===========================================================================
## Log transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
	  signature="logarithm",
	  function(expr,envir,enclos)
	  {    
            function(df)
              {   
                parameter=slot(expr,"parameters")
                if(class(parameter)!="transformReference")
                  {
                    parameter=eval(parameter)(df)
                  }
                else
                  {
                    parameter=resolveTransformReference(parameter,df)
                  }

                temp=slot(expr,"a")*parameter
                result=vector(mode="numeric",length=length(temp))
                result[temp>0]=log(temp[temp>0])*slot(expr,"b")
                return(result)
	      }
	  }
          )  

## ===========================================================================
## Exponential transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature="exponential",
          function(expr,envir,enclos)
          {    function(df)
                 {   
                   parameter=slot(expr,"parameters")
                   if(class(parameter)!="transformReference")
                     {
                       parameter=eval(parameter)(df)
                     }
                   else
                     {
                       parameter=resolveTransformReference(parameter,df)
                     }   
                   result=exp(parameter/slot(expr,"b"))/slot(expr,"a")
                 }
             }
          )

## ===========================================================================
## Inverse hyperbolic sin transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature="asinht",
          function(expr,envir,enclos)
          {    
            function(df)
              {
                parameter=slot(expr,"parameters")
                if(class(parameter)!="transformReference")
                  {
                    parameter=eval(parameter)(df)
                  }
                else
                  {
                    parameter=resolveTransformReference(parameter,df)
                  }   
                result=slot(expr,"b")*asinh(slot(expr,"a")*parameter)
              }
          }
          )

## ===========================================================================
## Hyperbolic sin transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature="sinht",
          function(expr,envir,enclos)
          {    function(df)
                 {
                   parameter=slot(expr,"parameters")
                   if(class(parameter)!="transformReference")
                     {
                       parameter=eval(parameter)(df)
                     }
                   else
                     {
                       parameter=resolveTransformReference(parameter,df)
                     }   
                   result=(sinh(parameter/
                                (slot(expr,"a"))))/(slot(expr,"b"))
                   return(result)
                 }
             }
          ) 

## ===========================================================================
## Hyperlog transformation 
## ---------------------------------------------------------------------------


setMethod("eval", 
	  signature="hyperlog",
	  function(expr,envir,enclos)
	  {    
            function(df)
              {  parameter=slot(expr,"parameters")
                 if(class(parameter)!="transformReference")
                   {
                     parameter=eval(parameter)(df)
                   }
                 else
                   {
                     parameter=resolveTransformReference(parameter,df)
                   }       
                 args=parameters
                 a=slot(expr,"a")
                 b=slot(expr,"b")
                 solveEH(args,a,b,df)
               }
	  }
          )

## ===========================================================================
## Splitscale transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature="splitscale",
          function(expr,envir,enclos)
          {    
            function(df)
              {   	
                parameter=slot(expr,"parameters")
                if(class(parameter)!="transformReference")
                  {
                    parameter=eval(parameter)(df)
                  }
                else
                  {
                    parameter=resolveTransformReference(parameter,df)
                  }     
                args=parameter
                
                transitionChannel=slot(expr,"transitionChannel")
                r=slot(object,"r")
                maxValue=slot(object,"maxValue")
                
                
                b=transitionChannel/2
                d=2*r*log10(2.71828)/transitionChannel
                log10t= -2*log10(2.71828)*r/transitionChannel +log10(maxValue)
                t=10^(log10t)
                a=transitionChannel/(2*t)
                
                log10ct= (a*t+b)*d/r
                c=(10^log10ct)/t
                
                idx <- which(args <= t)
                idx2 <- which(args > t)
                if(length(idx2)>0)
                  args[idx2] <- log10(c*args[idx2])*r/d
                if(length(idx)>0)
                  args[idx] <- a*args[idx]+b
                args
              }
	  }
          )
## ===========================================================================
## Inverse Splitscale transformation 
## ---------------------------------------------------------------------------
setMethod("eval", 
          signature="invsplitscale",
          function(expr,envir,enclos)
          {    
            function(df)
              {    
                parameter=slot(expr,"parameters")
                if(class(parameter)!="transformReference")
                  {
                    parameter=eval(parameter)(df)
                  }
                else
                  {
                    parameter=resolveTransformReference(parameter,df)
                  }        	
                                        #args=slot(object,"parameters")
                args=parameter
                transitionChannel=slot(expr,"transitionChannel")
                r=slot(expr,"r")
                maxValue=slot(expr,"maxValue")
                
                b=transitionChannel/2
                d=2*r*log10(2.71828)/transitionChannel
                log10t= -2*log10(2.71828)*r/transitionChannel +log10(maxValue)
                t=10^(log10t)
                a=transitionChannel/(2*t)
                log10ct= (a*t+b)*d/r
                c=(10^log10ct)/t            
                
                thresh=t*a+b
                args=parameter
                
                idx<-which(args<=thresh)
                idx2<-which(args>thresh)
                if(length(idx2)>0)
                  {
                    args[idx2]= (10^(args[idx2]*(d/r)))/c
                  }
                if(length(idx>0))
                  args[idx]=(args[idx]-b)/a
                args
                
              }
	  }
          )



## ===========================================================================
## Transformation reference
## ---------------------------------------------------------------------------
setMethod("eval",
          signature=signature(expr="transformReference",envir="missing",enclos="missing"),
	  function(expr,envir,enclos)
	  {
            (slot(expr,"searchEnv")[[slot(expr,"transformationId")]])
          }
	 )	
	
####----------------------------------------------------------------------------
#### 
#### Reference to a predefined gate 
#### 
####----------------------------------------------------------------------------

setMethod("eval",
          signature="filterReference",
	  function(expr,envir,enclos)
	  {   
              slot(expr,"env")[[slot(expr,"name")]]  #retrieved using name instead of the filterId
          }
	 )	
