## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## to flowSet
setMethod("[", c("flowSet"),
          function(x, i, j,..., drop=FALSE)
      {
          if(missing(i) && missing(j)) 
              return(x)
          if(!missing(j)){
              if(is.character(j))
                 colnames(x) <- colnames(x)[match(j, colnames(x))]
              else
                  colnames(x) <- colnames(x)[j] 
              if(any(is.na(colnames(x))))
                  stop("Subset out of bounds", call.=FALSE)
          }
          orig <- x@frames
          fr  <- new.env(hash=TRUE, parent=emptyenv())
          if(missing(i)) {
              for(nm in ls(orig))
                  fr[[nm]] <- orig[[nm]][, j, ..., drop=FALSE]
              pd <- phenoData(x)
          } else {
              if(is.numeric(i) || is.logical(i)) {
                  copy <- sampleNames(x)[i]
              } else {
                  copy <- i
                  i <- match(i,sampleNames(x))
              }
              if(any(is.na(copy)))
                  stop("Subset out of bounds", call.=FALSE)
              if(missing(j))
                  for(nm in copy)
                      fr[[nm]] <- orig[[nm]][, , ..., drop=FALSE]
              else
                  for(nm in copy)
                      fr[[nm]] <- orig[[nm]][, j, ..., drop=FALSE]
              pd <- phenoData(x)[i,]
          }
          fr <- as(fr,"flowSet")
          phenoData(fr) <- pd
          return(fr)
      })

## to flowFrame
setMethod("[[","flowSet",
          function(x, i, j, ...)
      {
          if(length(i) != 1)
              stop("subscript out of bounds (index must have length 1)")
          fr <- x@frames[[if(is.numeric(i)) sampleNames(x)[[i]] else i]]
          if(!missing(j))
              fr <- fr[,j]
          return(fr)
      })

## to flowFrame
setMethod("$", c("flowSet", "character"), function(x,name) x[[name]])



## ==========================================================================
## accessor and replace methods for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames","flowSet",function(x, do.NULL="missing",prefix="missing")
          x@colnames)


## fixme: should we make sure that parameter names are changed for each frame?
setReplaceMethod("colnames","flowSet",function(x,value) {
	x@colnames = value
	x
})



## ==========================================================================
## Allow for the extraction and replacement of phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("phenoData","flowSet",function(object) object@phenoData)
setMethod("phenoData<-","flowSet",function(object,value) {
	current <- phenoData(object)
	#Sanity checking
	if(nrow(current) != nrow(value))
		stop("phenoData must have the same number of rows as ",
                     "flow files")
	#Make sure all of the original frames appear in the new one.
	if(!all(sampleNames(current)%in%sampleNames(value)))
		stop("The sample names no longer match.")
	object@phenoData <- value
	object
})



## ==========================================================================
## directly access the pData data frame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("pData","flowSet",function(object) pData(object@phenoData))



## ==========================================================================
## set and extract the varLabels of the phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
## sampleNames method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("sampleNames", "flowSet", function(object) 
       sampleNames(phenoData(object)))

setReplaceMethod("sampleNames","flowSet",function(object,value)
             {
                 oldNames <- sampleNames(object)
                 if(length(oldNames)!=length(value) ||
                    !is.character(value))
                     stop(" replacement values must be character vector ",
                          "of length equal to number of frames in the set'")
                 env <- new.env(hash=TRUE,parent=emptyenv())
                 for(f in seq_along(oldNames))
                     assign(value[f], get(oldNames[f], object@frames), env)
                 pd <- phenoData(object)
                 sampleNames(pd) <- value
                 object@phenoData <- pd
                 object@frames <- env
                 return(object)
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
## accessor method for length of flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("length","flowSet",function(x) nrow(pData(phenoData(x))))



## ==========================================================================
## show method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show","flowSet",function(object) {
	cat("A flowSet with",length(object),"experiments.\n\n")
        if(any(varMetadata(phenoData(object))$labelDescription != "Name")){
            show(phenoData(object))
            cat("\n")
        }
	cat("  column names:\n  ")
	cat(paste(object@colnames,sep=","))
	cat("\n")
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
			phenoData(res) = phenoData(x)[sampleNames(x),]
		} else if(all(sapply(res,is.numeric)) && diff(range(sapply(res,length))) == 0) {
			res = do.call(rbind,res)
		}
	}
	res
})


## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
setMethod("compensate",signature("flowSet","matrix"),
          function(x,spillover,inv=TRUE, ...)
          fsApply(x,compensate,spillover,inv=inv, ...))



## ==========================================================================
## Transformation methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",signature=signature(`_data`="flowSet"),
          function(`_data`,...) {
              fsApply(`_data`,transform,...)
          })

setMethod("transform",signature(`_data`="missing"),function(...) {
    funs = list(...)
    io   = names(funs)
    ##Consistency check
    if(!all(sapply(funs,is.function)))
        stop("All transforms must be functions")
    if(!all(sapply(io,is.character)))
        stop("All transforms must be named")
    new("transformList",transforms=lapply(seq(along=funs),function(i)
                        new("transformMap",input=io[i],output=io[i],
                            f=funs[[i]])))
})



## ==========================================================================
## filter methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## for filters
setMethod("filter",signature=signature(x="flowSet", filter="filter"),
          function(x,filter)
      {
          if(!all(parameters(filter) %in% colnames(x)))
              stop("parameters in the filter definition don't ",
                   "match the parameters in the flowSet", call.=FALSE)
            fsApply(x,function(x) filter(x,filter))
      })

## for filterSets
setMethod("filter",signature(x="flowSet", filter="filterSet"),
          function(x, filter) fsApply(x, function(x) filter(x, filter)))

## for named lists of filters. Names of the list items have to correspond
## to sampleNames in the set.
setMethod("filter",signature(x="flowSet",filter="list"),
          function(x,filter)
      {
          if(is.null(names(filter)))
              stop("Filter list must have names to do something reasonable")
          nn <- names(filter)
          sn <- sampleNames(x)
          unused <- nn[!(nn %in% sn)]
          notfilter <-  sn[!(sn %in% nn)]
          ## Do some sanity checks
          if(length(unused) > 0)
              warning(paste("Some filters were not used:\n",
                            paste(unused, sep="", collapse=", ")),
                      call.=FALSE)
          if(length(notfilter) > 0)
              warning(paste("Some frames were not filtered:\n",
                            paste(notfilter, sep="", collapse=", ")),
                      call.=FALSE)
          common <- intersect(nn, sn)
          res <- list(length(common))
          for(f in common)
              res[[f]] <- filter(x[[f]], filter[[f]])
          return(res)
      })



## ==========================================================================
## split methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## split a flowSet by a single filter
setMethod("split",signature("flowSet","ANY"),
          function(x, f, drop=FALSE, population=NULL, prefix=NULL,
                   flowSet=FALSE, ...)
      {
          ## Split always returns a list
          sample.name <- sampleNames(x)
          fsApply(x,function(y) {
              l <- split(y, f, drop, population, prefix, flowSet=flowSet, ...)
              names(l) <- paste(names(l), "in", sample.name[1])
              sample.name <<- sample.name[-1]
              l
          }, simplify=FALSE)
      })

compatibleFilters <- function(f1, f2)
{
    if(class(f1) != class(f2))
        stop("Classes of filters don't match:\n", class(f1), " vs. ",
             class(f2), call.=FALSE)
    if(is(f1, "filterResult")){
        ff1 <- f1@filterDetails[[1]]$filter
        ff2 <- f2@filterDetails[[1]]$filter
        if(class(ff1) != class(ff2))
            stop("Classes of filters don't match:\n", class(ff1), " vs. ",
                 class(ff2), call.=FALSE)
        if(!(all(parameters(ff1) == parameters(ff2))))
            stop("Classes of filters don't match:\n",
                 paste(parameters(ff1), collapse=", "), " vs. ",
                 paste(parameters(ff2), collapse=", "), call.=FALSE)
    }
}

## split a flowSet according to a list of filters or filterResults
## of equal length
setMethod("split",signature("flowSet","list"),
          function(x, f, drop=FALSE, population=NULL,
                   prefix=NULL, flowSet=FALSE, ...)
      {
          ## A lot of sanity checking up front
          sample.name <- sampleNames(x)
          lf <- length(f)
          lx <- length(x)
          if(lf!=lx)
              stop("list of filterResults or filters must be same",
                   "length as flowSet")
          if(!all(sapply(f, is, "filter")))
              stop("Second argument must be list of filterResults or filters")
          lapply(f, compatibleFilters,  f[[1]])
          ## split everything or just some populations (if multipleFilterResult)
          if(is.null(population))
              if(!is.null(names(f[[1]])))
                  population <- names(f[[1]])
              else
                  population <- 1
          finalRes <- vector(mode="list", length=length(population))
          names(finalRes) <- population
          for(p in population){
              res <- vector(mode="list", length=lf)
              for(i in 1:lf){
                  l <- split(x[[i]], f[[i]], population=p,
                             prefix, flowSet=FALSE, ...)
                  res[[i]] <- l[[1]]
                  names(res)[i] <- paste(names(l), "in", sample.name[i])
              }
              if(flowSet)
                  finalRes[[p]] <- flowSet(res)
              else
                  finalRes[[p]] <- res
          }
          return(if(length(finalRes)==1) finalRes[[1]] else finalRes)
      })

## split a flowSet according to a factor, character or numeric 
setMethod("split",signature("flowSet","factor"),
          function(x,f,drop=FALSE,population=NULL,
                   prefix=NULL,flowSet=FALSE,...)
      {
          if(!is.atomic(f) || length(f)!=length(x))
              stop("split factor must be same length as flowSet") 
          gind <- split(1:length(f), f)
          res <- vector(mode="list", length=length(gind))
          for(g in seq_along(gind))
              res[[g]] <- x[gind[[g]]]
          return(res)
      })
setMethod("split",signature("flowSet","numeric"),
          function(x,f,drop=FALSE,population=NULL,
                   prefix=NULL,flowSet=FALSE,...)
          split(x, factor(f)))
setMethod("split",signature("flowSet","character"),
          function(x,f,drop=FALSE,population=NULL,
                   prefix=NULL,flowSet=FALSE,...)
          split(x, factor(f)))




## ==========================================================================
## Subset methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## by filter or filter result
setMethod("Subset", signature("flowSet","ANY"),
          function(x,subset,select,...)
      {
          y <- if(missing(select))
              fsApply(x, Subset, subset, ...)
          else
              fsApply(x, Subset, subset, select, ...)
          phenoData(y) <- phenoData(x)
          y
      })

setMethod("Subset", signature("flowSet", "list"),
          function(x, subset, select, ...)
      {
          if(is.null(names(subset)))
              stop("Filter list must have names to do something reasonable")
          nn <- names(subset)
          sn <- sampleNames(x)
          unused <- nn[!(nn %in% sn)]
          notfilter <- sn[!(sn %in% nn)]
          ##Do some sanity checks
          if(length(unused) > 0)
              warning(paste("Some filters were not used:\n",
                            paste(unused,sep="",collapse=", ")), call.=FALSE)
          if(length(notfilter) > 0)
              warning(paste("Some frames were not filtered:\n",
                            paste(notfilter,sep="",collapse=", ")),
                      .call=FALSE)	
          if(length(x) != length(subset))
              stop("You must supply a list of the same length as the flowSet.")
          used <- nn[nn %in% sn]
          res <- as(structure(if(missing(select))
                              lapply(used, function(i) Subset(x[[i]],
                                                              subset[[i]],...))
          else
                              lapply(used, function(i)
                                     Subset(x[[i]], subset[[i]], select, ...)),
                              names=sampleNames(x)), "flowSet")
          phenoData(res) <- phenoData(x)
          return(res)
      })




## ==========================================================================
## rbind method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("rbind2",signature("flowSet","missing"), function(x, y) x)
setMethod("rbind2",signature("flowSet", "flowSet"),
          function(x, y)
      {
          env <- new.env(hash=TRUE, parent=emptyenv())
          lx <- sampleNames(x)
          ly <- sampleNames(y)
          if(any(lx %in% ly))
              stop("These flowSets contain overlapping samples.")
          for(i in lx)
              assign(i, x[[i]], env=env)
          for(i in ly)
              assign(i, y[[i]], env=env)
          pd1 <- phenoData(x)
          pd2 <- phenoData(y)
          if(!all(varLabels(pd1) == varLabels(pd2)))
              stop("The phenoData of the two frames doesn't match.",
                   call.=FALSE)
          fs <- as(env,"flowSet")
          pData(pd1) <- rbind(pData(pd1), pData(pd2))
          phenoData(fs) <- pd1
          return(fs)
      })

setMethod("rbind2",signature("flowSet","flowFrame"),
          function(x,y)
      {
          ## create dummy phenoData
          pd <- phenoData(x)[1,]
          sampleNames(pd)
          pData(pd)[1,] <- NA
          tmp <- as(y, "flowSet")
          sampleNames(pd) <- sampleNames(tmp) <- "anonymous frame"
          phenoData(tmp) <- pd
          rbind2(x, tmp)
      })

setMethod("rbind2",signature("flowFrame","flowSet"),
          function(x,y) rbind2(y,x))
    





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


