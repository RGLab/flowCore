## ==========================================================================
## Allow for the extraction and replacement of phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("phenoData","flowSet",function(object) object@phenoData)
setMethod("phenoData<-","flowSet",function(object,value) {
	current = phenoData(object)
	#Sanity checking
	if(nrow(current) != nrow(value))
		stop("phenoData must have the same number of rows as flow files")
	#Make sure all of the original frames appear in the new one.
	if(!all(sampleNames(current)%in%sampleNames(value)))
		stop("The sample names no longer match.")
	object@phenoData = value
	object
})

setMethod("varLabels",
          signature=signature(object="flowSet"),
          function(object) varLabels(phenoData(object)))

setReplaceMethod("varLabels",
                 signature=signature(
                   object="flowSet",
                   value="ANY"),
                 function(object, value) {
                     pd <- phenoData(object)
                     varLabels(pd) <- value
                     object@phenoData <- pd
                     object
                 })

## ==========================================================================
## accessor and replace methods for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames","flowSet",function(x, do.NULL="missing",prefix="missing")
          x@colnames)

setReplaceMethod("colnames","flowSet",function(x,value) {
	x@colnames = value
	x
})


## ==========================================================================
## accessor method for length of flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("length","flowSet",function(x) nrow(pData(phenoData(x))))


## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setMethod("[",c("flowSet"),function(x,i,j,...,drop=FALSE) {
    if(missing(drop)) drop = FALSE
	if(missing(i) && missing(j)) 
		return(x)
	if(!missing(j))
		colnames(x) = colnames(x)[j]
	orig= x@frames
	fr  = new.env(hash=TRUE,parent=emptyenv())
	if(missing(i)) {
		for(nm in ls(orig)) fr[[nm]] = orig[[nm]][,j,...,drop=drop]
		pd = phenoData(x)
	} else {
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
		pd = phenoData(x)[i,]
	}
    	fr = as(fr,"flowSet")
	phenoData(fr) = pd
	fr
})
setMethod("[[","flowSet",function(x,i,j,...) {
	if(length(i)!=1)
		stop("subscript out of bounds (index must have length 1)")
	fr = x@frames[[if(is.numeric(i)) sampleNames(x)[[i]] else i]]
	fr
})

setMethod("$",c("flowSet","character"),function(x,name) x[[name]])


## ==========================================================================
## apply method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("fsApply",signature("flowSet","ANY"),function(x,FUN,...,simplify=TRUE,
                                                        use.exprs=FALSE) {
	if(missing(FUN))
		stop("fsApply function missing")
	FUN = match.fun(FUN)
	if(!is.function(FUN))
		stop("This is not a function!")
	## row.names and sampleNames had damn well better match, use this to
        ## give us access to the phenoData
	res = structure(lapply(sampleNames(x),function(n) {
		y = as(x[[n]],"flowFrame")
		FUN(if(use.exprs) exprs(y) else y,...)
	}),names=sampleNames(x))
	if(simplify) {
		if(all(sapply(res,is,"flowFrame"))) {
			res = as(res,"flowSet")
			phenoData(res) = phenoData(x)[sampleNames(x),]
		} else if(all(sapply(res,is.numeric)) && diff(range(sapply(res,length))) == 0) {
			res = do.call(rbind,res)
		}
	}
	res
})


## ==========================================================================
## Subset method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("Subset", signature("flowSet","ANY"),function(x,subset,select,...) {
	y = if(missing(select))
			fsApply(x,Subset,subset,...)
		else
			fsApply(x,Subset,subset,select,...)
	phenoData(y) = phenoData(x)
	y
})


setMethod("Subset", signature("flowSet","list"),function(x,subset,select,...) {
	if(is.null(names(subset)))
		stop("Filter list must have names to do something reasonable")
	nn = names(subset)
	sn = sampleNames(x)
	unused    = nn[!(nn %in% sn)]
	notfilter = sn[!(sn %in% nn)]
	#Do some sanity checks
	if(length(unused) > 0)
		warning(paste("Some filters were not used: ",paste(unused,sep=","),collapse=""))
	if(length(notfilter) > 0)
		warning(paste("Some frames were not filtered: ",paste(notfilter,sep=","),collapse=""))	
	if(length(x) != length(subset))
		stop("You must supply a list of the same length as the flowSet.")
	used = nn[nn %in% sn]
	res = as(structure(
		if(missing(select))
			lapply(used,function(i) Subset(x[[i]],subset[[i]],...))
		else
			lapply(used,function(i) Subset(x[[i]],subset[[i]],select,...)),
		names=sampleNames(x)),"flowSet")
	phenoData(res) = phenoData(x)
	res
})

## ==========================================================================
## split method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("split",signature("flowSet","ANY"),function(x,f,drop=FALSE,population=NULL,prefix=NULL,...) {
	#Split always returns a list
	sample.name = sampleNames(x)
	fsApply(x,function(y) {
		l = split(y,f,drop,population,prefix,...)
		names(l) = paste(names(l),"in",sample.name[1])
		sample.name <<- sample.name[-1]
		l
	},simplify=FALSE)
})


## ==========================================================================
## keyword method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("keyword",signature("flowSet","list"),function(object,keyword) {
	do.call("data.frame",c(lapply(keyword,function(k) {
		I(sapply(sampleNames(object),function(n) keyword(object[[n]],k)))
	}),list(row.names=sampleNames(object))))
})


## ==========================================================================
## rbind method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("rbind2",signature("flowSet","missing"),function(x,y) x)
setMethod("rbind2",signature("flowSet","flowSet"),function(x,y) {
	env = new.env(hash=TRUE,parent=emptyenv())
	lx  = sampleNames(x)
	ly  = sampleNames(y)
	if(any(lx %in% ly))
		stop("These flowSets contain overlapping samples.")
	for(i in lx) assign(i,x[[i]],env=env)
	for(i in ly) assign(i,y[[i]],env=env)
	fs            = as(env,"flowSet")
	pData(phenoData(fs)) = rbind(pData(phenoData(x)),pData(phenoData(y)))
	fs
})

setMethod("rbind2",signature("flowSet","flowFrame"),function(x,y) {
	##Should be able to find a flowFrame
})


## ==========================================================================
## show method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show","flowSet",function(object) {
	cat("A flowSet with",length(object),"experiments.\n\n")
	show(phenoData(object))
	cat("  \n  column names:\n  ")
	cat(paste(object@colnames,sep=","))
	cat("\n")
})


## ==========================================================================
## sampleNames method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("sampleNames", "flowSet", function(object) 
       sampleNames(phenoData(object)))


## ===========================================================================
## spillover method
## ---------------------------------------------------------------------------
setMethod("spillover","flowSet",function(x,unstained=NULL,patt=NULL,fsc="FSC-A",
                                         ssc="SSC-A",method="median") {
	if(is.null(unstained)) {
		stop("Sorry, we don't yet support unstained cells blended with stained cells")
	} else {
		## We often only want spillover for a subset of the columns 
		allcols = colnames(x)
		cols    = if(is.null(patt)) allcols else grep(patt,allcols,value=TRUE)
		
		
		if(is.numeric(fsc)) fsc = allcols[fsc]
		if(is.numeric(ssc)) ssc = allcols[ssc]

		

		if(is.na(match(fsc,allcols)))
			stop("Could not find forward scatter parameter. Please set the fsc parameter")
		if(is.na(match(ssc,allcols)))
			stop("Could not find side scatter parameter. Please set the ssc parameter")
		#Ignore these guys if they somehow got into cols.
		cols = cols[-match(c(fsc,ssc),cols)]

		## There has got to be a better way of doing this...
		if(!is.numeric(unstained)) {
			unstained = match(unstained,sampleNames(x))
			if(is.na(unstained))
				stop("Baseline not in this set.")
		}
		#Check to see if the unstained sample is in the list of stains. If not, we need
		#to add it, making it the first row and adjust the unstained index accordingly.
		#If it is there we adjust to the appropriate index.
		n2f   = norm2Filter(fsc,ssc,scale.factor=1.5)
		inten = fsApply(Subset(x,n2f),each_col,method)[,cols]
		inten = pmax(sweep(inten[-unstained,],2,inten[unstained,]),0)
		inten = sweep(inten,1,apply(inten,1,max),"/")
		row.names(inten) = colnames(inten)[apply(inten,1,which.max)]
		inten[colnames(inten),]
	}
})


## ===========================================================================
## compensate methods
## ---------------------------------------------------------------------------
setMethod("compensate",signature("flowSet","matrix"),
          function(x,spillover) fsApply(x,compensate,spillover))


## ==========================================================================
## Transformation methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",signature=signature(`_data`="flowSet"),function(`_data`,...) {
	fsApply(`_data`,transform,...)
})
setMethod("transform",signature(`_data`="missing"),function(...) {
	funs = list(...)
	io   = names(funs)
	#Consistency check
	if(!all(sapply(funs,is.function)))
		stop("All transforms must be functions")
	if(!all(sapply(io,is.character)))
		stop("All transforms must be named")
	new("transformList",transforms=lapply(seq(along=funs),function(i)
                              new("transformMap",input=io[i],output=io[i],f=funs[[i]])))
})


## ==========================================================================
## filter method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("filter",signature=signature(x="flowSet",filter="filter"),
          function(x,filter) {
            fsApply(x,function(x) filter(x,filter))
})

setMethod("filter",signature(x="flowSet",filter="list"),function(x,filter) {
	if(is.null(names(filter)))
		stop("Filter list must have names to do something reasonable")
	nn = names(filter)
	sn = sampleNames(x)
	unused    = nn[!(nn %in% sn)]
	notfilter = sn[!(sn %in% nn)]
	#Do some sanity checks
	if(length(unused) > 0)
		warning(paste("Some filters were not used: ",paste(unused,sep=","),collapse=""))
	if(length(notfilter) > 0)
		warning(paste("Some frames were not filtered: ",paste(notfilter,sep=","),collapse=""))
})
