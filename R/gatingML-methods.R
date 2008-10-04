resolveTransformReference<-function(trans,df)
{   
  if(class(trans)=="transformReference")
    {   trans=eval(trans)
        trans=resolveTransformReference(trans,df)
      }
  else
    {   
      trans=eval(trans)(df)
    }
}


EH<-function(y,a,b,argument)
{
  if(y>=0)
    return(10^(y/a)+b*y/a-1-argument)
  else
    return(-10^(-y/a)+b*y/a+1-argument)
}

solveEH<-function(args,a,b,df)
{
  temp=eval(args)(df)
  len=length(temp)
  result=0;
  while(len!=0)
    {
      result[len]=uniroot(EH,c(-100,100),tol=0.001,a=a,b=b,temp[len])$root
      len=len-1
    }
  return(result)
                                        #return(uniroot(EH,c(-10000,10000),tol=0.001,a=a,b=b,eval(args[[1]])(df)))
}


## ===========================================================================
## Evaluates all the transformations in the parameter slot of the filter 
## and returns a list containing the modified flowFrame and a modified filter 
## with slot parameter referencing the columns in the modified flowframe
## ---------------------------------------------------------------------------
resolveTransforms<-function(x,filter)
{

  data=exprs(x)
  recCoerce <- function(y){
    if(is(y, "filterReference")){
      y <- as(y, "concreteFilter")
      recCoerce(y)
    }else
    y
  }  
  filter <- recCoerce(filter)
  if(class(filter)!="intersectFilter" & class(filter)!="complementFilter" &
     class(filter)!="unionFilter" & class(filter)!="expressionFilter" &
     class(filter)!="subsetFilter")
    {
      parameters <- filter@parameters
      len <- length(parameters)
      charParam <- list()
      
      while(len>0)
        {
          if(class(parameters[[len]])!="unitytransform") 
            {   
              ## process all transformed parameters
              charParam[[len]]=sprintf("_NEWCOL%03d_",len) 

              if (class(eval(parameters[[len]]))=="compensatedParameter") 
                {   ##deals with compensated parameter
                  newCol=eval(eval(parameters[[len]]))(x)
                }
              else if(class(parameters[[len]])=="transformReference")    
                {   ##deals with transform references
                  temp=parameters[[len]]
                  newCol=resolveTransformReference(parameters[[len]],data)
                }
              else                           
                {   ##deals with all other transforms
                  newCol=eval(parameters[[len]])(data)
                }
              
              newCol=as.matrix(newCol)
              colnames(newCol)=sprintf("_NEWCOL%03d_",len) 
              data <- cbind(data, newCol)
            } 
          else
            {
              charParam[[len]]=slot(parameters[[len]],"parameters")
            }                
          len=len-1                               
        }
      
      parameters(filter) <- charParam
    } 
  
  y <- flowFrame(data)
  identifier(y) <- identifier(x)
  return(list(y,filter))
}
