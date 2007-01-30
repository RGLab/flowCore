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
                                                        data=data.frame(name=I(frameList),row.names=frameList),varMetadata=data.frame(labelDescription="Name",row.names="name")))
},function(from,value) {
    if(!canCoerce(value,"AnnotatedDataFrame"))
      stop("Must be able to coerce 'value' to an AnnotatedDataFrame for use as metadata for this set")
    from            = as(from,"flowSet")
    phenoData(from) = as(value,"AnnotatedDataFrame")
    from
})
## ==========================================================================
                                        

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
## ==========================================================================


## ==========================================================================
## Allow for the extraction and replacement of phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("phenoData","flowSet",function(object) object@phenoData)
setMethod("phenoData<-","flowSet",function(object,value) {
	current = phenoData(object)
	#Sanity checking
	if(nrow(current) != nrow(value))
		stop("phenoData must have the same number of rows as flow files")
	#If the row.names have changed (not just reordered we should remap them)
	if(!all(sampleNames(current)==sampleNames(value)))
		stop("The sample names no longer match.")
	object@phenoData = value
	object
})
## ==========================================================================


## ==========================================================================
## accessor method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames","flowSet",function(x, do.NULL="missing",prefix="missing")
          x@colnames)
## ==========================================================================

## ==========================================================================
## replace method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setReplaceMethod("colnames","flowSet",function(x,value) {
	x@colnames = value
	x
})
## ==========================================================================


## ==========================================================================
## accessor method for length of flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("length","flowSet",function(x) nrow(pData(phenoData(x))))
## ==========================================================================


## ==========================================================================
## subsetting method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setMethod("[",c("flowSet"),function(x,i,j,...,drop=FALSE) {
	if(missing(drop)) drop = FALSE
	if(missing(i) && missing(j)) 
		return(x)
	if(!missing(j))
		colnames(x) = colnames(x)[j]
	orig= x@frames
	fr  = new.env(hash=TRUE,parent=emptyenv())
	if(missing(i))
		for(nm in ls(orig)) fr[[nm]] = orig[[nm]][,j,...,drop=drop]
	else {
		if(is.numeric(i) || is.logical(i)) {
			copy = phenoData(x)$name[i]
		} else {
			copy = i
			i    = match(i,phenoData(x)$name)
		}
		if(missing(j))
			for(nm in copy) fr[[nm]] = orig[[nm]][,,...,drop=drop]
		else
			for(nm in copy) fr[[nm]] = orig[[nm]][,j,...,drop=drop]
	}
	new("flowSet",frames = fr,phenoData=phenoData(x)[i,],colnames=if(missing(j)) x@colnames else x@colnames[j])
})
setMethod("[[","flowSet",function(x,i,j,...) {
	if(length(i)!=1)
		stop("subscript out of bounds (index must have length 1)")
	fr = x@frames[[if(is.numeric(i)) sampleNames(x)[[i]] else i]]
	fr
})

## ==========================================================================


## ==========================================================================
## apply method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("fsApply",signature("flowSet","ANY"),function(x,FUN,...,simplify=TRUE,use.exprs=FALSE) {
	FUN = match.fun(FUN)
	if(!is.function(FUN))
		stop("This is not a function!")
	# row.names and sampleNames had damn well better match, use this to give us access to the phenoData
	res = structure(lapply(sampleNames(x),function(n) {
		y = as(x[[n]],"flowFrame")
		FUN(if(use.exprs) exprs(y) else y,...)
	}),names=sampleNames(x))
	if(simplify) {
		if(all(sapply(res,is,"flowFrame"))) {
			res = as(res,"flowSet")
			if(all(sampleNames(res) == sampleNames(x)))
				phenoData(res) = phenoData(x)
		} else if(all(sapply(res,is.numeric)) && diff(range(sapply(res,length))) == 0) {
			res = do.call(rbind,res)
		}
	}
	res
})

setMethod("Subset", signature("flowSet","ANY"),function(x,subset,select,...) {
	if(missing(select))
		fsApply(x,Subset,subset,...)
	else
		fsApply(x,Subset,subset,select,...)
})

setMethod("split",signature("flowSet","ANY"),function(x,f,drop=FALSE,population=NULL,prefix=NULL,...) {
	#Split always returns a list
	fsApply(x,function(y) {
		l = split(y,f,drop,population,prefix,...)
		names(l) = paste(names(l),"in",sample.name)
		l
	},simplify=FALSE)
})

setMethod("keyword",signature("flowSet","list"),function(object,keyword) {
	do.call("data.frame",c(lapply(keyword,function(k) {
		I(sapply(sampleNames(object),function(n) keyword(object[[n]],k)))
	}),list(row.names=sampleNames(object))))
})

setMethod("rbind2",signature("flowSet","missing"),function(x,y) x)
setMethod("rbind2",signature("flowSet","flowSet"),function(x,y) {
	env = new.env(hash=TRUE,parent=emptyenv())
	lx  = sampleNames(x)
	ly  = sampleNames(y)
	if(lx %in% ly)
		stop("These flowSets contain overlapping samples.")
	for(i in lx) assign(i,x[[i]],env=env)
	for(i in ly) assign(i,y[[i]],env=env)
	fs            = as(env,"flowSet")
	phenoData(fs) = rbind2(phenoData(x),phenoData(y))
	fs
})


## ==========================================================================
## show method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show","flowSet",function(object) {
	cat("A flowSet with ",length(object)," experiments.\n\n")
	show(phenoData(object))
	cat("\nColumn names:\n")
	cat(paste(object@colnames,sep=","))
	cat("\n")
})
## ==========================================================================

setMethod("sampleNames", "flowSet", function(object) 
       sampleNames(phenoData(object)))
