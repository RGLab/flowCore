## last modified Nov 3, 2006
## pdh: the method exprs wasn't being passed by the package to the workspace. I'm not sure why. So I
##      changed to use the slot call directly

## ==========================================================================
## filtering function for polygonal gates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rectangleFiltering <- function(filter,flowObject,parent){
    parameters=filter@parameters
    min=filter@min
    max=filter@max
    ## test for validity of objects
    
    if(class(flowObject)!="flowFrame")
      stop("flowObject must be of class flowFrame.")
    if (is.null(flowObject@exprs)) 
        stop("There is no data to createFilter.")
    
    data=flowObject@exprs
    
    
    ## ncells is the size of the original data
    ncells = nrow(data)
    ## a 1 means the cell (row) was in the filter, a 0 means it wasn't
    selectNew = rep(0,ncells)
    ## check for the parent. If it is a filterResult get the parentSet
    parentSet = getParentSet(parent,ncells)
    ## a 1 means the cell (row) was in the filter, a 0 means it wasn't
    selectCurrent = (parentSet == 1)
    data = data[selectCurrent,parameters,drop=FALSE]
    ## npop is the size of the population actually used
    npop = nrow(data)
    
    rangeSel = rep(TRUE,npop)
    ndim = length(parameters)

    ## note that the data matrix has already been reordered to correspond with the
    ## dimensions specified by the gate
    for(i in 1:ndim){
      rangeSel <- rangeSel & ((data[,i] >= min[i]) & (data[,i] <  max[i]))
    }

    newIndex <-  ((1:ncells)[selectCurrent])[rangeSel]
    selectNew[newIndex] <- 1
    return(selectNew)
}

## ==========================================================================
## filtering function for polygonal gates
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
polygonFiltering <- function(filter,flowObject,parent){
    parameters=filter@parameters
    boundaries=filter@boundaries
    ## test for validity of objects


    if(class(flowObject)!="flowFrame")
      stop("flowObject must be of class flowFrame.")
    if (is.null(flowObject@exprs)) 
      stop("There is no data to createFilter.")
    
    data=flowObject@exprs

  
    ## ncells is the size of the original data
    ncells = nrow(data)
    ## a 1 means the cell (row) was in the filter, a 0 means it wasn't
    selectNew = rep(0,ncells)
    ## check for the parent. If it is a filterResult get the parentSet
    parentSet = getParentSet(parent,ncells)
    ## a 1 means the cell (row) was in the filter, a 0 means it wasn't
    selectCurrent = (parentSet == 1)
    data = data[selectCurrent,parameters,drop=FALSE]
    ## npop is the size of the population actually used
    npop = nrow(data)
    
    rangeSel = rep(TRUE,npop)
    ndim = length(parameters)
    
    
    if (ndim==1) {
        rangeSel <- rangeSel & ((data[,parameters] >= boundaries[1,]) & (data[,parameters] <  boundaries[2,]))
    }
    if (ndim==2) {
        ##lib.there <- require(prada)
        ##if (lib.there==FALSE){
        ##    stop("Library prada can not be found but is needed")
        ##}            
        data <- cbind(as.numeric(data[,parameters[1]]),as.numeric(data[,parameters[2]]))
        rangeSel <-as.logical(.Call("inPolygon",data, boundaries,PACKAGE="flowCore"))
    }
    newIndex <-  ((1:ncells)[selectCurrent])[rangeSel]
    selectNew[newIndex] <- 1
    return(selectNew)
}

## ==========================================================================
## filtering function with bivariate normal distribution 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
normFiltering <- function(filter,flowObject,parent){
    parameters=filter@parameters
    scale.factor=filter@scale.factor
    method=filter@method
    ## test for validity of objects


    if(class(flowObject)!="flowFrame")
      stop("flowset must be of class flowFrame.")
    if (is.null(flowObject@exprs)) 
      stop("There is no data to createFilter.")
    
    data=flowObject@exprs
    
    
    if(length(parameters)!=2) {
        stop("Two data dimensions must be specified.")
    }
    ## ncells is the size of the original data
    ncells = nrow(data)
    ## a 1 means the cell (row) was in the filter, a 0 means it wasn't
    selectNew = rep(0,ncells)
    ## check for the parent. If it is a filterResult get the parentSet
    parentSet = getParentSet(parent,ncells)
    ## a 1 means the cell (row) was in the filter, a 0 means it wasn't
    selectCurrent = (parentSet == 1)
    data = data[selectCurrent,parameters,drop=FALSE]
    ## npop is the size of the population actually used
    npop = nrow(data)
    
    rangeSel = rep(TRUE,npop)
    ndim = length(parameters)
    
    
    ## lib.there <- require(prada)
    ##if (lib.there==FALSE){
    ##    stop("Library prada can not be found but is needed")
    ##}
    fn.fit <- fitNorm2(data[,parameters[1]],data[,parameters[2]],scalefac=scale.factor,method=method)
       
    
    newIndex <-  ((1:ncells)[selectCurrent])[fn.fit$sel]
    selectNew[newIndex] <- 1
    ans <- list(sel=selectNew,mu=fn.fit$mu,S=fn.fit$S)
    return(ans)
}

## ==========================================================================
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getParentSet <- function(parent,ncells) {
 	if(missing(parent)) {
   		## if the parent is missing assume that it is all cells
		parentSet = rep(1,ncells)
 	}
	else {  	
		if(class(parent) != "filterResult") {
			stop("parent must be of class filterResult.")
		}
		else {
		if(is.null(parent@subSet)) {
				stop("The filterResult parent must have a subSet slot.")
			}
			else {
				parentSet = parent@subSet
			}  	
		}
	}
	parentSet
}
 
