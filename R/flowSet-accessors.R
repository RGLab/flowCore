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
	new("flowSet",frames = fr,phenoData=phenoData(x)[i,],colnames=if(missing(j))
            x@colnames else x@colnames[j])
})
setMethod("[[","flowSet",function(x,i,j,...) {
	if(length(i)!=1)
		stop("subscript out of bounds (index must have length 1)")
	fr = x@frames[[if(is.numeric(i)) sampleNames(x)[[i]] else i]]
	fr
})


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
			if(all(sampleNames(res) == sampleNames(x)))
				phenoData(res) = phenoData(x)
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
	if(missing(select))
		fsApply(x,Subset,subset,...)
	else
		fsApply(x,Subset,subset,select,...)
})


## ==========================================================================
## split method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("split",signature("flowSet","ANY"),function(x,f,drop=FALSE,population=NULL,prefix=NULL,...) {
	#Split always returns a list
	fsApply(x,function(y) {
		l = split(y,f,drop,population,prefix,...)
		names(l) = paste(names(l),"in",sample.name)
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
		"is.null(unstained)"
	} else {
		## We often only want spillover for a subset of the columns 
		cols = if(is.null(patt)) colnames(x) else grep(patt,colnames(x)
                  ,value=TRUE)
		## Ignore the forward and sidescatter channels if they managed to get included
		if(!is.na(match(fsc,cols))) cols = cols[-match(fsc,cols)] 
		if(!is.na(match(ssc,cols))) cols = cols[-match(ssc,cols)]
		## There has got to be a better way of doing this...
		stains = phenoData(x)$name
		stains = stains[-match(unstained,stains)]
		
		## Grab the baseline from the unstained values
		baseline = apply(exprs(x[[unstained]])[x[[unstained]] %in%
                  norm2Filter(fsc,ssc,scale.factor=1.5),cols],2,method)
		
		## Now do the same thing to all the stains and sweep out the baseline
                ## to figure out the motion on all of the channels.
		inten = sweep(sapply(stains,function(s)
			apply(exprs(x[[s]])[x[[s]] %in%
                                            norm2Filter(fsc,ssc,scale.factor=1.5),cols],
                              2,method)),2,baseline)
		## We assume that the highest intensity channel in each column is the signal
                ## channel. If something weird happens here, you probably screwed up your comp
                ## controls (or you're using an awful channel like PacO in which case the mean
                ## is probably the recommended statistic).
		colnames(inten) = cols[apply(inten,2,which.max)]
		#Now normalize row-wise to figure out the % spillover and ensure that
                ## any negative values are set to 0 since negative compensation is even more
                ## insane than negative fluoresence.
		inten = apply(inten,2,function(x) ifelse(x<0,0,x/max(x)))		
		t(inten[,match(cols,colnames(inten))])
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
