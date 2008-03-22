## ==========================================================================
## We can convert a factor, logical or a numeric into a filterResult by
## selecting a specific filterResult type. This is done through the standard
## R coercion techniques.
## --------------------------------------------------------------------------
setAs("factor", "filterResult",function(from)
      new("multipleFilterResult",filterId="",subSet=from))

setAs("logical","filterResult",function(from)
      new("logicalFilterResult",filterId="",subSet=from))

setAs("numeric","filterResult",function(from)
      new("randomFilterResult",filterId="",subSet=from))



## ==========================================================================
## We can also convert some filterResult types directly to logical types,
## though in general it is not possible. We provide the means for logical 
## and random filter types. 
## --------------------------------------------------------------------------
setAs("filterResult","logical",function(from)
      stop("Unable to convert to a logical vector"))

setAs("logicalFilterResult","logical",function(from)
      from@subSet)

setAs("randomFilterResult","logical",function(from)
      runif(length(from@subSet))<from@subSet)



## ==========================================================================
## Allow the coercion of resolvable filters (i.e. those derived from
## filterResult) to be composed and then converted into a logical vector.
## This allows for a lot of processing to be done simply using the filter
## results.
## --------------------------------------------------------------------------
setAs("filter","logical",function(from)
      stop("Only resolved filters can be converted to a logical vector."))

setAs("subsetFilter","logical",function(from)
      as(from@filters[[1]],"logical") & as(from@filters[[2]],"logical"))
setAs("intersectFilter","logical",function(from)
      
      apply(sapply(from@filters,as,Class="logical"),1,all))
setAs("unionFilter","logical",function(from)
      apply(sapply(from@filters,as,Class="logical"),1,any))

setAs("complementFilter","logical",function(from)
      !as(from@filters[[1]],"logical"))



## ==========================================================================
## Allows for the resolution of filterReferences and formulas
## --------------------------------------------------------------------------
setAs("filterReference","concreteFilter",function(from) {
    x = from@env[[from@name]]
    if(is.null(x)) stop(paste("Unable to resolve filter reference:",
                              from@name))
    x
})

setAs("formula","filter",function(from) {
    f = as(from[[length(from)]],"filter")
    if(length(from) == 3 && from[[2]] != ".")
        f@filterId = as.character(from[[2]])
    f
})

setAs("character","filter",function(from) {
    filterReference(as(find(from,mode='S4'),"environment"),from)
})

setAs("name","filter",function(from) as(as.character(from),"filter"))
setAs("call","filter",function(from) {
    filters = lapply(from[-1],as,Class="filter")
    eval(as.call(c(from[[1]],filters)))
})



## ==========================================================================
## These exist primarily to support making copies of filterSet objects
## --------------------------------------------------------------------------
setAs("filterReference","call",function(from) as.symbol(from@name))

setAs("filter","call",function(from) {
    nam  = names(getSlots(class(from)))
    vals = structure(lapply(nam,function(n) {
        v = slot(object=from,n)
        if(is.call(v)) as.call(c(as.symbol("quote"),v)) else v
    }),names=nam)
    as.call(c(as.symbol("new"),class(from),vals))
})

setAs("subsetFilter","call",function(from) {
    eval(as.call(c(as.symbol('~'),as.symbol(from@filterId),
                   as.call(c(as.symbol("%subset%"),
                             lapply(from@filters,as,Class="call"))))))
})


## Helper function for converting lists into binary call trees
binaryHelper = function(op,l) {
	x = l[[1]]
	op = as.symbol(op)
	for(i in 2:length(l)) x = as.call(c(op,x,l[[i]]))
	x
}

setAs("intersectFilter","call",function(from) {
    eval(as.call(c(as.symbol('~'),as.symbol(from@filterId),
                   binaryHelper('&',lapply(from@filters,as,Class="call")))))
})

setAs("unionFilter","call",function(from) {
    eval(as.call(c(as.symbol('~'),as.symbol(from@filterId),
                   binaryHelper('|',lapply(from@filters,as,Class="call")))))
})

setAs("complementFilter","call",function(from) {
    if(length(from@filters) > 1)
        stop("Whoops. Complements only work on one filter right now.")
    eval(as.call(c(as.symbol('~'),as.symbol(from@filterId),
                   as.call(c(as.symbol("!"),as(from@filters[[1]],"call"))))))
})



## ==========================================================================
## Make ourselves some copies!
## --------------------------------------------------------------------------
setAs("filterSet","list",function(from)	{
    nam = ls(env=from@env)
    out = lapply(nam,function(n) as(from[[n]],"call"))
    ##   actually want names everywhere.
    ##	nam[sapply(out,is,"formula")] = ""
    names(out) = nam
    out
})

setAs("list","filterSet",function(from) {
    fs = filterSet()
    e  = new.env()
    n    = names(from)
    if(is.null(n) || length(n)==0) n = rep("",length(from))
    for(i in 1:length(from)) {
        filter = from[[i]]
        name   = if(!is.null(n[i]) && nzchar(n[i])) n[i] else NULL
        fs[[name]] = eval(filter,e)
    }
    fs
})

## ==========================================================================
## Coerce method: Convert an environment to a flowSet.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("environment","flowSet",function(from) {
    frameList  = ls(env=from)
    isFrame    = sapply(frameList,function(f) is(get(f,env=from),"flowFrame"))
    if(!all(isFrame))
      warning("Some symbols are not flowFrames.",
              "They will be ignored but left intact.")
    ## If specified, remove extraneous symbols from the environment
    ## before continuing
    frameList = frameList[isFrame]
    
    ##Check the column names
    colNames = sapply(frameList,function(f) colnames(from[[f]]))
    if(is.null(dim(colNames)))
        dim(colNames) <- c(ncol(from[[frameList[[1]]]]), length(frameList))
    if(!all(apply(colNames,2,"==",colNames[,1])))
        stop("Column names for all frames do not match.")
    new("flowSet",frames=from,colnames = colNames[,1],
        phenoData=new("AnnotatedDataFrame",
        data=data.frame(name=I(frameList),row.names=frameList),
        varMetadata=data.frame(labelDescription="Name",row.names="name")))
})

## ==========================================================================
## Convert a list to a flowSet by creating an environment and coerce THAT
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("list","flowSet",function(from) {
    if(is.null(names(from)))
        names(from) = paste("V",seq(1,length(from)),sep="")
    as(l2e(from,new.env(hash=T,parent=emptyenv())),"flowSet")
})

## ==========================================================================
## Convert a flowSet to a list
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("flowSet","list",function(from) {
    ret <- list()
    for(i in sampleNames(from))
        ret[[i]] <- from[[i]]
    return(ret)
})


## ==========================================================================
## Coerce a flowFrame to a flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("flowFrame", "flowSet", function(from) {
    flowSet(from)
})




## ==========================================================================
## Coerce a flowSet to a flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("flowSet","flowFrame",function(from) {
    if(length(from) == 1)
        from[[1]]
    else {
        ## The parameters need to match all frames
        params = parameters(from[[1]])
        allParams <- fsApply(from, function(x) pData(parameters(x))$name)
        if(!all(sapply(allParams, identical, pData(params)$name)))
            stop("parameters must be the same for all frames")
        ## making sure we are not doing too many copies of the data
        lens <- fsApply(from, nrow)
        exp <- matrix(ncol=nrow(params)+1, nrow=sum(lens))
        colnames(exp) <- c(colnames(from), "Original")
        offset <- 1
        for(i in 1:length(from)){
            rows <- offset:(offset+lens[i,]-1)
            exp[rows, 1:nrow(params)] <- exprs(from[[i]])
            exp[rows,"Original"] <- rep(i, lens[i,])
            offset <- offset+lens[i,]
        }
        pData(params)["Original",c("name", "desc")] <- c("Original",
                              "Original frame")
        desc  <- list(description="Synthetic Frame",
                      sampleNames=sampleNames(from))
        new("flowFrame",exprs=exp,parameters=params,description=desc)
    }
})





setAs("filterSummary","data.frame",function(from) {
	data.frame("true"=from@true,"false"=from@count-from@true,"count"=from@count,"p"=from@p,"q"=1-from@q,row.names=from@name)
})




