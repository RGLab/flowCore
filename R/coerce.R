## ==========================================================================
## We can convert a factor, logical or a numeric into a filterResult by
## selecting a specific filterResult type. This is done through the standard
## R coercion techniques.
## --------------------------------------------------------------------------
setAs("factor", "filterResult",function(from)
      new("multipleFilterResult",parameters=character(0),filterId="",subSet=from))
setAs("logical","filterResult",function(from)
      new("logicalFilterResult",parameters=character(0),filterId="",subSet=from))
setAs("numeric","filterResult",function(from)
      new("randomFilterResult",parameters=character(0),filterId="",subSet=from))

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
## Coerce method: Convert an environment to a flowSet.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("environment","flowSet",function(from) {
    frameList  = ls(env=from)
    isFrame    = sapply(frameList,function(f) is(get(f,env=from),"flowFrame"))
    if(!all(isFrame))
      warning("Some symbols are not flowFrames. They will be ignored but left intact.")
    ##If specified, remove extraneous symbols from the environment before continuing
    frameList = frameList[isFrame]
    
    ##Check the column names
    colNames = sapply(frameList,function(f) colnames(from[[f]]))
    if(!all(apply(colNames,2,"==",colNames[,1])))
      stop("Column names for all frames do not match.")
    new("flowSet",frames=from,colnames = colNames[,1],phenoData=new("AnnotatedDataFrame",
        data=data.frame(name=I(frameList),row.names=frameList),
        varMetadata=data.frame(labelDescription="Name",row.names="name")))
  },function(from,value) {
    if(!canCoerce(value,"AnnotatedDataFrame"))
      stop("Must be able to coerce 'value' to an AnnotatedDataFrame for use as metadata for this set")
    from            = as(from,"flowSet")
    phenoData(from) = as(value,"AnnotatedDataFrame")
    from
  })
                                        

## ==========================================================================
## Coerce method: Convert a list to a flowSet by creating an environment and converting THAT
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("list","flowSet",function(from) {
    env = new.env(hash=TRUE,parent=emptyenv())
    multiassign(from,env=env)
	as(env,"flowSet")
    },function(from,value) {
	env = new.env(hash=TRUE,parent=emptyenv())
	multiassign(from,env=env)
	as(env,"flowSet") <- value
})
