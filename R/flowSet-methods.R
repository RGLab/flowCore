## ==========================================================================
## Convert an environment to a flowSet.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("environment","flowSet",function(from) {
	frameList  = ls(env=from)
	isFrame    = sapply(frameList,function(f) is(get(f,env=from),"flowFrame"))
	if(!all(isFrame))
		warning("Some symbols are not flowFrames. They will be ignored but left intact.")
	#If specified, remove extraneous symbols from the environment before continuing
	frameList = frameList[isFrame]

	#Check the column names
	colNames = sapply(frameList,function(f) colnames(get(f,env=from)))
	if(!all(apply(colNames,2,"==",colNames[,1])))
		stop("Column names for all frames do not match.")
	phenoData = new("AnnotatedDataFrame",
		data=data.frame(name=I(frameList)),varMetadata=data.frame(labelDescription="Name of FCS file",row.names="name"))
	for(f in frameList){
		x = get(f,env=from)
		colnames(x) = NULL
		assign(f,x,env=from)
	}
	new("flowSet",frames=from,phenoData=phenoData,colnames=colNames[,1])
})
## ==========================================================================
                                        

## ==========================================================================
## Convert a list to a flowSet by creating an environment and converting THAT
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs("list","flowSet",function(from) {
	env = new.env(hash=TRUE,parent=emptyenv())
	multiassign(from,env=env)
	as(env,"flowSet")
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
	if(length(unique(value$name)) != nrow(value))
		stop("phenoData must have a name column to uniquely identify flowFrames.")
	#If the names are different, we need to remap the environment assuming that the ordering is the same.
	if(!all(current$name==value$name)) {
		fr = object@frames
		for(i in seq(along=current$name)) {
			#If the name differs then swap the names
			if(current$name[i] != value$name[i]) {
				x = get(current$name[i],env=fr)
				rm(current$name[i],env=fr)
				assign(x,value$name[i],env=fr)
			}
		}
	}
	object@phenoData = value
	object
})
## ==========================================================================


## ==========================================================================
## accessor method for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames","flowSet",function(x, do.NULL="missing",prefix="missing") x@colnames)
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
	fr = x@frames[[if(is.numeric(i)) phenoData(x)$name[i] else i]]
	colnames(exprs(fr)) = x@colnames
	fr
})

## ==========================================================================

setMethod("fsApply",signature("flowSet","function"),function(x,FUN,...) {
	res = as(structure(lapply(phenoData(x)$name,function(n) FUN(x[[n]],...)),names=phenoData(x)$name),"flowSet")
	phenoData(res) = phenoData(x)
	res
})


## ==========================================================================
## show method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show","flowSet",function(object) {
	cat("A flowSet with ",length(object)," experiments.\n\n")
	cat("Column names:\n")
	cat(paste(object@colnames,sep=","))
	cat("\n")
})
## ==========================================================================

setMethod("sampleNames", "flowSet", function(object) 
       sampleNames(phenoData(object))
