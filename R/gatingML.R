## recursively resolve transformation references until we get to a
## realized transform object
resolveTransformReference<-function(trans,df)
{  
    if(is(trans=="transformReference")){
        trans <- eval(trans,trans@searchEnv)
        trans <- resolveTransformReference(trans, df)
    }else{  
        if(!is(trans, "function")) 
            trans <- eval(trans)(df)
    }
    trans
}


## FIXME NG: Please document
EH <- function(y, a, b, argument)
{
    if(y>=0)
        return(10^(y/a)+b*y/a-1-argument)
    else
        return(-10^(-y/a)+b*y/a+1-argument)
}

solveEH<-function(args,a,b,df)
{
  temp <- eval(args)(df)
  len <- length(temp)
  result <- 0
  while(len!=0)
  {
      result[len]=uniroot(EH,c(-100,100),tol=0.001,a=a,b=b,temp[len])$root
      len <- len-1
  }
  return(result)
}


## ===========================================================================
## Evaluates all the transformations in the parameter slot of the filter 
## and returns a list containing the modified flowFrame and a modified filter 
## with slot parameter referencing the columns in the modified flowframe
## ---------------------------------------------------------------------------
resolveTransforms<-function(x,filter)
{    
    ## Filters need to be fully realized for this to work and we have to
    ## resolve all filterReferences
    recCoerce <- function(y){
        if(is(y, "filterReference")){
            y <- as(y, "concreteFilter")
            recCoerce(y)
        }else
        y
    }  
    filter <- recCoerce(filter)
    if(is(filter, "parameterFilter")){
        parameters <- filter@parameters
        charParam <- list()
        for(len in seq_along(parameters)){
            if(!is(parameters[[len]], "unitytransform")){   
                ## process all transformed parameters
                charParam[[len]] <- sprintf("_NEWCOL%03d_", len) 
                if (is(eval(parameters[[len]]), "compensatedParameter")){
                    ##deals with compensated parameter
                    newCol <- eval(eval(parameters[[len]]))(x)
                }else if(is(parameters[[len]], "transformReference")){
                    ##deals with transform references
                    temp <- parameters[[len]]
                    newCol <- resolveTransformReference(parameters[[len]], data)
                }else{
                    ##deals with all other transforms
                    newCol <- eval(parameters[[len]])(data)
                }               
                newCol <- as.matrix(newCol)
                colnames(newCol) <- sprintf("_NEWCOL%03d_", len) 
                x <- cbind(x, newCol)
            }else
               charParam[[len]] <- parameters(parameters[[len]])
        }
        parameters(filter) <- charParam
    } 
    return(list(data=x, filter=filter))
}
