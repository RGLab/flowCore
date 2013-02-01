## recursively resolve transformation references until we get to a
## realized transform object
resolveTransformReference<-function(trans,df)
{  
    if(is(trans, "transformReference")){
        trans <- eval(trans)
        if(!is(trans, "compensatedParameter"))
          trans <- resolveTransformReference(trans, df)
        else
          trans <- eval(trans)(df)
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

solveEH<-function(args,a,b)
{  
    #temp=eval(args)(df)
    len=length(args)
    result=0;
    
    while(len!=0)
    {
        step<-2
        checkme<-EH(-step,a,b,args[len])*EH(step,a,b,args[len])
        while(checkme>0)
        {
          step<-step*2
          checkme<-EH(-step,a,b,args[len])*EH(step,a,b,args[len])
        } 
        result[len]=uniroot(EH,c(-step,step),tol=0.001,a=a,b=b,args[len])$root
        len=len-1
            
    }

    return(matrix(result))
}

## ===========================================================================
## Evaluates all the transformations in the parameter slot of the filter 
## and returns a list containing the modified flowFrame and a modified filter 
## with slot parameter referencing the columns in the modified flowframe
## ---------------------------------------------------------------------------
resolveTransforms <- function(x, filter)
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
                if (is(parameters[[len]], "compensatedParameter")){
                    ##deals with compensated parameter
                    ##newCol <- eval(eval(parameters[[len]]))(x)
                    newCol <- resolveTransformReference(parameters[[len]],
                                                        exprs(x))
                }else if(is(parameters[[len]], "transformReference")){
                    ##deals with transform references
                    newCol <- resolveTransformReference(parameters[[len]],
                                                        exprs(x))
                }else{
                    ##deals with all other transforms
                    newCol <- eval(parameters[[len]])(exprs(x))
                }               
                newCol <- as.matrix(newCol)
                colnames(newCol) <- sprintf("_NEWCOL%03d_", len) 
                x <- cbind2(x, newCol)
            }else
               charParam[[len]] <- parameters(parameters[[len]])
        }
        parameters(filter) <- charParam
    } 
    return(list(data=x, filter=filter))
}
