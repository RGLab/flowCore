## ==========================================================================
## flowSets are basically lists flowFrames
## ==========================================================================






## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## to flowSet
#' @export
setMethod("[",
		signature=signature(x="flowSet"),
		definition=function(x, i, j, ..., drop=FALSE)
		{
			if(missing(i) && missing(j)) 
				return(x)
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
			if(!missing(j)){
				if(is.character(j))
					colnames(fr) <- colnames(x)[match(j, colnames(x))]
				else
					colnames(fr) <- colnames(x)[j] 
				if(any(is.na(colnames(fr))))
					stop("Subset out of bounds", call.=FALSE)
			}
			return(fr)
		})

## to flowFrame
#' @export
setMethod("[[",
		signature=signature(x="flowSet"),
		definition=function(x, i, j, ...)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have length 1)")
			fr <- x@frames[[if(is.numeric(i)) sampleNames(x)[[i]] else i]]
			if(!missing(j))
				fr <- fr[,j]
			return(fr)
		})

## to flowFrame
#' @export
setMethod("$",
          signature=signature(x="flowSet"),
          definition=function(x, name) x[[name]])

## replace a flowFrame
setReplaceMethod("[[",
		signature=signature(x="flowSet",
				value="flowFrame"),
		definition=function(x, i, j, ..., value)
		{
			if(length(i) != 1)
				stop("subscript out of bounds (index must have ",
						"length 1)")
                        cnx <- colnames(x)
                        cnv <- colnames(value)
                        if(length(cnx) != length(cnv) || !all(sort(cnv) == sort(cnx)))
                            stop("The colnames of this flowFrame don't match ",
                                 "the colnames of the flowSet.")
                        
			sel <- if(is.numeric(i)) sampleNames(x)[[i]] else i
			x@frames[[sel]] <- value
			return(x)
		})

## ==========================================================================
## accessor and replace methods for slot colnames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @noRd
#' @export
setMethod("colnames",
		signature=signature(x="flowSet"),
		definition=function(x, do.NULL="missing", prefix="missing")
			x@colnames)

#' @export
setReplaceMethod("colnames",
		signature=signature(x="flowSet",
				value="ANY"),
		definition=function(x, value)
		{
			x@colnames <- value
                        for(i in sampleNames(x))
                            colnames(x@frames[[i]]) <- value
			x
		})

#' @rdname markernames
#' @export
setMethod("markernames",
          signature=signature(object="flowSet"),
          definition=function(object){
            
          res <- lapply(object@frames, function(fr){
              markernames(fr)
            })
          
          res <- unique(res)
          if(length(res) > 1)
            warning("marker names are not consistent across samples within flowSet")
          else
            res <- res[[1]]
          res
        })
#' @rdname markernames
#' @export
setReplaceMethod("markernames",
                 signature=signature(object="flowSet", value="ANY"), function(object, value){
                   for(i in ls(object@frames))
                     markernames(object@frames[[i]]) <- value
                   object
         })
## ==========================================================================
## Allow for the extraction and replacement of phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("phenoData",
		signature=signature(object="flowSet"),
		definition=function(object) object@phenoData)

#' @export
setMethod("phenoData<-",
		signature=signature(object="flowSet",
				value="ANY"),
		definition=function(object, value)
		{
			current <- phenoData(object)
			## Sanity checking
			if(nrow(current) != nrow(value))
				stop("phenoData must have the same number of rows as ",
						"flow files")
			## Make sure all of the original frames appear in the new one.
			if(!all(sampleNames(current)%in%sampleNames(value)))
				stop("The sample names no longer match.")
            #validity check for 'name' column
            df <- pData(value)
            if(!"name" %in% colnames(df))
              pData(value)[["name"]] = rownames(df)
            
              
			object@phenoData <- value
			object
		})



## ==========================================================================
## directly access the pData data frame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("pData",
		signature=signature(object="flowSet"),
		definition=function(object) pData(phenoData(object)))

#' @export
setReplaceMethod("pData",
		signature=signature(object="flowSet",
				value="data.frame"),
		definition=function(object,value)
		{
			pd <- phenoData(object)
			pData(pd) <- value
			phenoData(object) <- pd
			object
		})

## ==========================================================================
## set and extract the varLabels of the phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("varLabels",
		signature=signature(object="flowSet"),
		function(object) varLabels(phenoData(object)))

#' @export
setReplaceMethod("varLabels",
		signature=signature(object="flowSet",
				value="ANY"),
		definition=function(object, value)
		{
			pd <- phenoData(object)
			varLabels(pd) <- value
			object@phenoData <- pd
			object
		})

#' @export
setMethod("varMetadata",
		signature=signature(object="flowSet"),
		definition=function(object) varMetadata(phenoData(object)))

#' @export
setReplaceMethod("varMetadata",
		signature=signature(object="flowSet",
				value="ANY"),
		definition=function(object, value)
		{
			pd <- phenoData(object)
			varMetadata(pd) <- value
			object@phenoData <- pd
			object
		})



## ==========================================================================
## sampleNames method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("sampleNames",
		signature=signature(object="flowSet"),
		definition=function(object) 
			sampleNames(phenoData(object)))

## Note that the replacement method also replaces the GUID for each flowFrame
#' @export
setReplaceMethod("sampleNames",
		signature=signature(object="flowSet"),
		definition=function(object, value)
		{
			oldNames <- sampleNames(object)
			value <- as.character(value)
			if(length(oldNames)!=length(value) ||
					!is.character(value))
				stop(" replacement values must be character vector ",
						"of length equal to number of frames in the set'",
						call.=FALSE)
			if(any(duplicated(value)))
				stop("Replacement values are not unique.", call.=FALSE)
			env <- new.env(hash=TRUE,parent=emptyenv())
			for(f in seq_along(oldNames)){
				tmp <- get(oldNames[f], object@frames)
				identifier(tmp) <- value[f]
				assign(value[f], tmp, env)
			}
			pd <- phenoData(object)
			sampleNames(pd) <- value
			object@phenoData <- pd
			object@frames <- env
			return(object)
		})


## ==========================================================================
## keyword method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("keyword",
		signature=signature(object="flowSet",
				keyword="list"),
		definition=function(object, keyword)
		{
			keys <-  fsApply(object, function(x) unlist(keyword(x, keyword)))
			if(!is.null(dim(keys))){
				colnames(keys) <- gsub("\\..*$", "", colnames(keys))
				rownames(keys) <- sampleNames(object)
			}
			return(keys)
		})

#' @export
setMethod("keyword",
		signature=signature(object="flowSet",
				keyword="ANY"),
		definition=function(object, keyword)
			keyword(object, as.list(keyword)))

setReplaceMethod("keyword", signature=c("flowSet", "list"),
		definition=function(object, value){
			for(i in seq_along(value)){
				vals <- rep(value[[i]], length(object))
				for(j in seq_len(length(object))){
					thisVal <- list(vals[[j]])
					names(thisVal) <- names(value)[i] 
					keyword(object[[j]]) <- thisVal
				}
			}
			object
		})



## ==========================================================================
## apply method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' Apply a Function over values in a flowSet
#' 
#' \code{fsApply}, like many of the \code{apply}-style functions in R, acts as an
#' iterator for \code{flowSet} objects, allowing the application of a function
#' to either the \code{flowFrame} or the data matrix itself. The output can then
#' be reconstructed as either a \code{flowSet}, a list, or a matrix depending on
#' options and the type of objects returned.
#' 
#' @name fsApply
#' @aliases fsApply fsApply,flowSet,ANY
#' 
#' @usage 
#' fsApply(x, FUN, \dots, simplify=TRUE, use.exprs=FALSE)
#' 
#' @param x \code{\link[flowCore:flowSet-class]{flowSet}} to be used
#' @param FUN the function to be applied to each element of \code{x}
#' @param simplify logical (default: TRUE); if all true and all objects are
#' \code{flowFrame} objects, a \code{flowSet} object will be constructed. If
#' all of the values are of the same type there will be an attempt to construct
#' a vector or matrix of the appropriate type (e.g. all numeric results will
#' return a matrix).
#' @param use.exprs logical (default: FALSE); should the \code{FUN} be applied
#' on the \code{\link[flowCore:flowFrame-class]{flowFrame}} object or the
#' expression values.
#' @param \dots optional arguments to \code{FUN}.
#' 
#' @author B. Ellis
#' @seealso \code{\link{apply}}, \code{\link{sapply}}
#' @keywords iteration
#' @examples
#' 
#' fcs.loc <- system.file("extdata",package="flowCore")
#' file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")
#' samp <- read.flowSet(file.location[1:3]) 
#' 
#' #Get summary information about each sample.
#' fsApply(samp,summary)
#' 
#' #Obtain the median of each parameter in each frame.
#' fsApply(samp,each_col,median)
#' 
#' 
#' @export
setMethod("fsApply",
		signature=signature(x="flowSet",
				FUN="ANY"),
		definition=function(x,FUN,...,simplify=TRUE, use.exprs=FALSE)
		{
			if(missing(FUN))
				stop("fsApply function missing")
			FUN <- match.fun(FUN)
			if(!is.function(FUN))
				stop("This is not a function!")
			## row.names and sampleNames had damn well better match, use this to
			## give us access to the phenoData
			res <- structure(lapply(sampleNames(x),function(n) {
								# y <- as(x[[n]],"flowFrame")
			          y <- x[[n]]
								FUN(if(use.exprs) exprs(y) else y,...)
							}),names=sampleNames(x))
			if(simplify) {
				if(all(sapply(res,is,"flowFrame"))) {
					res <- as(res,"flowSet")
					phenoData(res) = phenoData(x)[sampleNames(x),]
				} else if(all(sapply(res,is.numeric)) || all(sapply(res,is.character)) &&
						diff(range(sapply(res,length))) == 0) {
					res <- do.call(rbind,res)
				}
			}
			res
		})


## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
#' @export
setMethod("compensate",
		signature=signature(x="flowSet",
				spillover="ANY"),
		definition=function(x, spillover)
			fsApply(x, compensate, spillover))

#' @export
setMethod("compensate",
          signature=signature(x="flowSet",
              spillover="data.frame"),
          definition=function(x, spillover)
            selectMethod("compensate"
                    , signature=signature(x="flowSet",spillover="ANY"))(x, spillover)
      )
#' @export
setMethod("compensate",
    signature=signature(x="flowSet",
        spillover="list"),
    definition=function(x, spillover){
      samples <- sampleNames(x)
      if(!all(samples %in% names(spillover)))
        stop("names of the compensation list must match the sample names of 'flowSet'!")
      
      res <- structure(lapply(samples, function(sn)compensate(x[[sn]], spillover[[sn]])),names=sampleNames(x))
      res <- as(res,"flowSet")
      phenoData(res) = phenoData(x)[sampleNames(x),]
      res
    })
      

## ==========================================================================
## Transformation methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("transform",
		signature=signature(`_data`="flowSet"),
		definition=function(`_data`, translist, ...)
		{
		  if(missing(translist))
		    fsApply(`_data`,transform, ...)
		  else if(is(translist, "transformList"))
			  fsApply(`_data`,transform, translist = translist, ...)
		  else if(is(translist, "list")){
		    sns <- sampleNames(`_data`)
		    if(!setequal(sns, names(translist)))
		      stop("names of 'translist' must be consistent with flow data!")
		    fs <- copyFlowSet(`_data`)
		    for(sn in sns)
		      fs[[sn]] <- transform(fs[[sn]], translist[[sn]])
		    fs
		  }else
		    stop("expect the second argument as a 'transformList' object or a list of 'transformList' objects!")
		  
		})

#' @export
setMethod("transform",
		signature=signature(`_data`="missing"),
		definition=function(...)
		{
			funs <- list(...)
			io <- names(funs)
			## Consistency check
			if(!all(sapply(funs,is.function)))
				stop("All transforms must be functions")
			if(!all(sapply(io,is.character)))
				stop("All transforms must be named")
			new("transformList",
					transforms=lapply(seq(along=funs),function(i)
                                                          new("transformMap",input=io[i],
                                                              output=io[i],
                                                              f=funs[[i]])))
		})



## ==========================================================================
## filter methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## These methods apply single filters or lists of filters to a
## flowSet object. In all cases, the output of the filtering operation is
## a filterResultList
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## for filters
#' @noRd
#' @export
setMethod("filter",
		signature=signature(x="flowSet",
				filter="filter"),
		definition=function(x, filter, method = "missing", sides = "missing", circular = "missing", init = "missing")
		{
			if(!all(parameters(filter) %in% colnames(x)))
				stop("parameters in the filter definition don't ",
						"match the parameters in the flowSet", call.=FALSE)
			res <- fsApply(x,function(x) filter(x,filter))
			return(new("filterResultList", .Data=res, frameId=sampleNames(x),
							filterId=identifier(filter)))
		})

#' @noRd
#' @export
setMethod("filter",
		signature=signature(x="flowSet",
				filter="list"),
		definition=function(x, filter, method = "missing", sides = "missing", circular = "missing", init = "missing")
      {
          filter(x, filterList(filter))
      })

## for named lists of filters. Names of the list items have to correspond
## to sampleNames in the set. Filters in the filter list that can't be
## matched are ignored, for those that are missing, an "empty" dummy
## filterResult is produced
#' @noRd
#' @export
setMethod("filter",
		signature=signature(x="flowSet",
				filter="filterList"),
		definition=function(x, filter, method = "missing", sides = "missing", circular = "missing", init = "missing")
      {
          if(is.null(names(filter)))
              stop("'filter' must be a named list, where names correspond",
                   " to sample names in the flowSet", call.=FALSE)
          nn <- names(filter)
          sn <- sampleNames(x)
          unused <- nn[!(nn %in% sn)]
          notfilter <-  setdiff(sn, nn)
          ## Check for non-matching filters
          if(length(unused) > 0)
              warning(paste("Some filters were not used:\n",
                            paste(unused, sep="", collapse=", ")),
                      call.=FALSE)
          common <- intersect(nn, sn)
          res <- vector("list", length(x))
          fid <- character(length(x))
          names(res) <- names(fid) <- sampleNames(x)
          ## use all matching filters first
          for(f in common){
              res[[f]] <- filter(x[[f]], filter[[f]])
              fid[f] <- identifier(filter[[f]])
          }
          ## use dummy filters for all the rest (if any)
          if(length(notfilter)){
              warning(paste("Some frames were not filtered:\n",
                            paste(notfilter, sep="", collapse=", ")),
                      call.=FALSE)
              exp <- paste("rep(length(", parameters(x[[1]], names=TRUE)[1],
                           "))", sep="")
              dummyFilter <- char2ExpressionFilter(exp, filterId="dummy")
              res[notfilter] <- filter(x[notfilter], dummyFilter)
              fid[notfilter] <- identifier(dummyFilter)
          }
          return(new("filterResultList", .Data=res, frameId=sampleNames(x),
                     filterId=fid))
      })





## ==========================================================================
## Subset methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## by filter or filter result
#' @export
setMethod("Subset",
		signature=signature(x="flowSet",
				subset="ANY"),
		definition=function(x,subset,select,...)
		{
			y <- if(missing(select))
						fsApply(x, Subset, subset, ...)
					else
						fsApply(x, Subset, subset, select, ...)
			phenoData(y) <- phenoData(x)
			y
		})

#' @export
setMethod("Subset",
		signature=signature(x="flowSet",
				subset="filterResultList"),
		definition=function(x, subset, select, ...)
		{
			flowCore:::validFilterResultList(subset, x, strict=FALSE)
			res <- as(structure(if(missing(select))
										lapply(names(subset), function(i) Subset(x[[i]],
															subset[[i]],...))
									else
										lapply(names(subset), function(i)
													Subset(x[[i]], subset[[i]], select, ...)),
							names=sampleNames(x)), "flowSet")
			phenoData(res) <- phenoData(x)
			return(res)
		})


#' @export
setMethod("Subset",
		signature=signature(x="flowSet",
				subset="list"),
		definition=function(x, subset, select, ...)
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
#' @export
setMethod("rbind2",
		signature=signature(x="flowSet",
				y="missing"),
		definition=function(x, y) x)

#' @export
setMethod("rbind2",
		signature=signature(x="flowSet",
				y="flowSet"),
		definition=function(x, y)
		{
			env <- new.env(hash=TRUE, parent=emptyenv())
			lx <- sampleNames(x)
			ly <- sampleNames(y)
			if(any(lx %in% ly))
				stop("These flowSets contain overlapping samples.")
			for(i in lx)
				assign(i, x[[i]], envir = env)
			for(i in ly)
				assign(i, y[[i]], envir = env)
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

#' @export
setMethod("rbind2",
		signature=signature(x="flowSet",
				y="flowFrame"),
		definition=function(x,y)
		{
			## create dummy phenoData
			pd <- phenoData(x)[1,]
			
			pData(pd)[1,] <- NA
			tmp <- as(y, "flowSet")
            pData(pd)[1, "name"] <- sampleNames(pd) <- sampleNames(tmp) <- "anonymous frame"
			phenoData(tmp) <- pd
			rbind2(x, tmp)
		})

#' @export
setMethod("rbind2",
		signature=signature(x="flowFrame",
				y="flowSet"),
		definition=function(x,y) rbind2(y,x))






## ===========================================================================
## spillover method
## ---------------------------------------------------------------------------
#' Compute a spillover matrix from a flowSet
#' 
#' Spillover information for a particular experiment is often obtained by
#' running several tubes of beads or cells stained with a single color that can
#' then be used to determine a spillover matrix for use with
#' \code{\link{compensate}}.\cr\cr
#' When matching stain channels in \code{x} with the compensation controls, we
#' provide a few options. If \code{ordered}, we assume the ordering of the
#' channels in the flowSet object is the same as the ordering of the
#' compensation-control samples. If \code{regexpr}, we use a regular expression
#' to match the channel names with the names of each of the compensation control
#' \code{flowFrame}s (that is, \code{sampleNames(x)}, which will typically be the 
#' filenames passed to \code{\link{read.FCS}}).
#' By default, we must "guess" based on the largest statistic for the
#' compensation control (i.e., the row).\cr\cr
#' Additionally, matching of channels to compensation control files can
#' be accomplished using the \code{\link{spillover_match}} method, which allows
#' the matches to be specified using a csv file. The \link{flowSet} returned
#' by the \code{spillover_match} method should then be used as the \code{x} argument
#' to \code{spillover} with \code{prematched = TRUE}.
#' 
#' 
#' The algorithm used is fairly simple. First, using the scatter parameters, we
#' restrict ourselves to the most closely clustered population to reduce the
#' amount of debris. The selected statistic is then calculated on all
#' appropriate parameters and the unstained values swept out of the matrix.
#' Every sample is then normalized to [0,1] with respect to the maximum value
#' of the sample, giving the spillover in terms of a proportion of the primary
#' channel intensity.
#' @name spillover-flowSet
#' @aliases spillover spillover,flowSet-method
#' 
#' @usage 
#' \S4method{spillover}{flowSet}(x, unstained = NULL, fsc = "FSC-A", 
#'                               ssc = "SSC-A", patt = NULL, method = "median", 
#'                               stain_match = c("intensity", "ordered", "regexpr"),
#'                               useNormFilt=FALSE, prematched = FALSE)
#' 
#' @param x A flowSet of compensation beads or cells
#' @param unstained The name or index of the unstained negative control
#' @param patt An optional regular expression defining which parameters should
#' be considered
#' @param fsc The name or index of the forward scatter parameter
#' @param ssc The name or index of the side scatter parameter
#' @param method The statistic to use for calculation. Traditionally, this has
#' been the median so it is the default. The mean is sometimes more stable.
#' @param stain_match Determines how the stain channels are matched with the
#' compensation controls. See details.
#' @param useNormFilt logical Indicating whether to apply a
#' \code{\link{norm2Filter}} first before computing the spillover
#' @param prematched a convenience argument specifying if the channels
#' have already been matched by spillover_match. This will override the
#' values of unstained and stain_match with unstained = "unstained" and
#' stain_match = "regexpr".
#' 
#' @return A matrix for each of the parameters
#' @author B. Ellis, J. Wagner
#' @seealso \code{\link{compensate}}, \code{\link{spillover_match}}
#' @references C. B. Bagwell & E. G. Adams (1993). Fluorescence spectral
#' overlap compensation for any number of flow cytometry parameters. in: Annals
#' of the New York Academy of Sciences, 677:167-184.
#' @keywords methods
#' @export
setMethod("spillover",
          signature = signature(x = "flowSet"),
          definition = function(x, unstained = NULL, fsc = "FSC-A",
                                ssc = "SSC-A", patt = NULL, method = "median",
                                stain_match = c("intensity", "ordered", "regexpr"),
                                useNormFilt = FALSE, prematched = FALSE) {
            if(prematched){
              unstained = "unstained"
              stain_match <- "regexpr"
            }
            stain_match <- match.arg(stain_match)
            
            if (is.null(unstained)) {
              stop("Sorry, we don't yet support unstained cells blended ",
                   "with stained cells", " please specify the name or index of unstained sample", call. = FALSE)
            } else {
              # Check validity of fsc, ssc, and convert to numeric
              if(!(is.numeric(fsc))){
                if(!(fsc %in% colnames(x))){
                  stop(paste0("Forward scatter parameter '", fsc, "' not found in this flowSet."), call. = FALSE)
                }else{
                  fsc <- match(fsc, colnames(x))
                }
              }else if(fsc > length(colnames(x)) || fsc < 1){
                stop(paste0("Forward scatter parameter index", fsc, "' is out of bounds."), call. = FALSE)
              }
              
              if(!(is.numeric(ssc))){
                if(!(ssc %in% colnames(x))){
                  stop(paste0("Side scatter parameter '", ssc, "' not found in this flowSet."), call. = FALSE)
                }else{
                  ssc <- match(ssc, colnames(x))
                }
              }else if(ssc > length(colnames(x)) || ssc < 1){
                stop(paste0("Side scatter parameter index", ssc, "' is out of bounds."), call. = FALSE)
              }
              
              ## We often only want spillover for a subset of the columns
              allcols <- colnames(x)
              cols <- if (is.null(patt)) {
                allcols
              } else {
                grep(patt, allcols, value = TRUE)
              }

              ## Ignore these guys if they somehow got into cols.
              cols <- cols[-c(fsc,ssc)]
              
              ## There has got to be a better way of doing this...
              if (!is.numeric(unstained)) {
                unstained <- match(unstained, sampleNames(x))
                if (is.na(unstained)) {
                  stop('Baseline(unstained sample) not found in this set. Note that the default
                        sample name searched for is "unstained". If this is incorrect, please 
                        specify the name or index of unstained sample.', call. = FALSE)
                }
              }

              ## Check to see if the unstained sample is in the list of
              ## stains. If not, we need to add it, making it the first
              ## row and adjust the unstained index accordingly.
              ## If it is there we adjust to the appropriate index.
              
              if (useNormFilt) {
                if (is.numeric(fsc)) {
                  fsc <- allcols[fsc]
                }
                if (is.numeric(ssc)) {
                  ssc <- allcols[ssc]
                }

                if (is.na(match(fsc, allcols))) {
                  stop("Could not find forward scatter parameter. ",
                       "Please set the fsc parameter", call. = FALSE)
                }
                if (is.na(match(ssc, allcols))) {
                  stop("Could not find side scatter parameter. ",
                       "Please set the ssc parameter", call. = FALSE)
                  n2f <- norm2Filter(fsc, ssc, scale.factor = 1.5)
                  x <- Subset(x, n2f)
                }
              }
              
              # Here, we match the stain channels with the compensation controls
              # if the user has specified it. Otherwise, we must "guess" below
              # based on the largest statistic for the compensation control
              # (i.e., the row).
              # If "ordered," we assume the ordering of the channels in the
              # flowSet object is the same as the ordering of the
              # compensation-control samples.
              # Another option is to use a regular expression to match the
              # channel names with the sampleNames of the flowSet.
              if (stain_match == "intensity") {
                channel_order <- NA
              } else if (stain_match == "ordered") {
                channel_order <- seq_along(sampleNames(x))[-unstained]
              } else if (stain_match == "regexpr") {
                channel_order <- sapply(cols, grep, x = sampleNames(x), fixed = TRUE)
                # Clip out those channels that do not match to a name
                matched <- channel_order[sapply(channel_order, length) != 0]
                if(anyDuplicated(matched)) {
                  stop("Multiple stains match to a common compensation-control name",
                       call. = FALSE)
                }
              }
              if(length(x) - 1 != length(cols))
              {
                stop("the number of single stained samples provided in this set doesn't match to the number of stained channels!")
              }
              if (method == "mode") {
                inten <- fsApply(x, function(flow_frame) {
                  modes <- sapply(cols, function(stain) {
                    density_stain <- density(exprs(flow_frame)[, stain])
                    with(density_stain, x[which.max(y)])
                  }, USE.NAMES = TRUE)
                  modes
                })
              } else {
                inten <- fsApply(x, each_col, method)[, cols]
              }

              # background correction
              inten <- pmax(sweep(inten[-unstained, ], 2, inten[unstained, ]), 0)

              # normalize by max of each row
              inten <- sweep(inten, 1, apply(inten, 1, max), "/")

              # Updates the "rownames" of the intensity matrix. If the channel
              # order was not set above, then a guess is made based on the
              # largest statistic for the compensation control (i.e., the row).
              if (any(is.na(channel_order))) {
                channel_order <- apply(inten, 1, which.max)
                if (anyDuplicated(channel_order) > 0) {
                  stop("Unable to match stains with controls based on intensity: ",
                       "a single stain matches to several multiple controls. ",
                       call. = FALSE)
                }
              }else{
                # Just bump-down the channel_order to account for
                # the removal of the unstained row. 
                channel_order <- channel_order - (channel_order > unstained)
                # Then reverse any shuffling so the rows are in the same
                # order as the channels for symmetry
                inten <- inten[order(channel_order),]
              }
              rownames(inten) <- colnames(inten)
              inten
          }
      })

## =========================================================================
## spillover_match method: Match channel names to compensation control filenames
## -------------------------------------------------------------------------
#' Construct a \code{flowSet} for use with \code{spillover} by matching channel
#' names to compensation control filenames
#' 
#' Spillover information for a particular experiment is often obtained by
#' running several tubes of beads or cells stained with a single color that can
#' then be used to determine a spillover matrix for use with
#' \code{\link{compensate}}.\cr\cr
#' This method facilitates construction of a \code{flowSet} of compensation
#' control \code{flowFrame}s using a simple file linking filenames to channels.
#' This resulting \code{flowSet} can then be used with \code{\link{spillover}}
#' using the option \code{prematched = TRUE}.\cr\cr
#' Matching stain channels to compensation controls is done via a csv
#' file (\code{matchfile}) with columns 'filename' and 'channel'. The 'channel' entries should
#' exactly match the channel names in the FCS files. The 'filename' should be
#' the FCS file name of each compensation control which should also be the
#' corresponding sample name in the \code{flowSet}. There should also be one unstained
#' control with the 'channel' entry of 'unstained'.\cr\cr
#' The method also allows for \code{x} to be missing if \code{path} is provided,
#' pointing to a directory containing the control FCS files.\cr\cr
#' 
#' 
#' @name spillover_match-flowSet
#' @aliases spillover_match spillover_match,flowSet-method
#' @usage 
#' \S4method{spillover_match}{flowSet}(x, fsc = "FSC-A", ssc = "SSC-A",
#'                                     matchfile = NULL, path)
#' @param x A flowSet of compensation beads or cells
#' @param fsc The name or index of the forward scatter parameter
#' @param ssc The name or index of the side scatter parameter
#' @param matchfile The name or path of the csv file holding the compensation control file
#' to channel matching information.
#' @param path The name or path of the directory containing the control FCS files to
#' be matched to channels by matchfile.
#' 
#' @return A \code{flowSet} with the sample names of its \code{flowFrames}
#' corresponding to the channels specified by the matchfile.
#' @author B. Ellis, J. Wagner
#' @seealso \code{\link{compensate}}, \code{\link{spillover}}
#' @keywords methods

#' @export
setMethod("spillover_match",
			 signature = signature(x = "flowSet"),
			 definition = function(x, fsc = "FSC-A", ssc = "SSC-A", matchfile = NULL, path) {
			 	if(!file.exists(matchfile)){
			 		stop("File ", matchfile, " not found. \n 
			 		     You must provide a file that maps single-stained 
               controls to channels (including one unstained control). \n 
			 		     It should have columns 'filename' and 'channel'", call. = FALSE)
			 	}
			   # Check validity of fsc, ssc, and convert to numeric
			  if(!(is.numeric(fsc))){
			    if(!(fsc %in% colnames(x))){
			      stop(paste0("Forward scatter parameter '", fsc, "' not found in this flowSet."), call. = FALSE)
			    }else{
			    fsc <- match(fsc, colnames(x))
			    }
			  }else if(fsc > length(colnames(x)) || fsc < 1){
			    stop(paste0("Forward scatter parameter index", fsc, "' is out of bounds."), call. = FALSE)
			  }
			   
			  if(!(is.numeric(ssc))){
			    if(!(ssc %in% colnames(x))){
			      stop(paste0("Side scatter parameter '", ssc, "' not found in this flowSet."), call. = FALSE)
			    }else{
			      ssc <- match(ssc, colnames(x))
			    }
			  }else if(ssc > length(colnames(x)) || ssc < 1){
			    stop(paste0("Side scatter parameter index", ssc, "' is out of bounds."), call. = FALSE)
			  }
			   
			 	match <- read.csv(matchfile, colClasses = "character")
			 	if(!all(sort(tolower(colnames(match)))[c(2,1)] == c("filename","channel"))){
			 		stop("File ", matchfile, 
               " should map single-stained controls to channels (including one unstained control). \n
			 		     It should have columns 'filename' and 'channel'",call. = FALSE)
			 	}
			 	if (sum(tolower(match$channel)%in%"unstained")==0) {
			 		stop("Sorry, we don't yet support unstained cells blended ",
			 			  "with stained cells", " please specify the name or index of unstained sample", call. = FALSE)
			 	} else if (sum(tolower(match$channel)%in%"unstained")>1) {
			 		stop("You may only specify one unstained sample in the match file.")
			 	} else {
			 	  
			 		unstained.match.idx <- which(tolower(match$channel)%in%"unstained")
			 		unstained <- as.character(match[unstained.match.idx,"filename"])
			 		
			 		allcols <- colnames(x)
			 		#check if the columns passed in via the match file agree with the columns in the fcs files.
			 		channels <- as.character(match$channel)
			 		if(!all(channels[-unstained.match.idx]%in%allcols)){
			 			stop("Not all channels in ", match, 
			 			     " match to channel names in the FCS files. Check your spelling and capitalization.")
			 		}
			 		if(sum(match$filename%in%unstained)!=1)
			 		  stop("Unique baseline (unstained sample) not found in this set.", call. = FALSE)
			 		
			 		# We often only want spillover for a subset of the columns
			 		# Drop those columns that don't have a line in the matchfile, while maintaining original column order
			 		keep <- c(fsc, ssc, match(channels[-unstained.match.idx], allcols))
			 		x <- x[,sort(keep)]
			 		
			 		# Make sure the sample names of x are given their corresponding channel names
			 		# (This is the guts of the file <-> channel matching logic)
			 		sampleNames(x) <- sapply(sampleNames(x), 
			 		                         function(to_match) 
			 		                         as.character(match[match$filename == to_match, "channel"]))
			 		x
			 	}
			 })

#' @export
#' @rdname spillover_match-flowSet
setMethod("spillover_match",
          signature = signature(x = "missing"),
          definition = function(x, fsc = "FSC-A", ssc = "SSC-A", matchfile, path){
          if(missing(path)){
            stop("If no flowSet is provided, a path must be provided to the control FCS files", call. = FALSE)
          }
          match <- read.csv(matchfile, colClasses = "character")  
          # By default, this will give the file names to the flowSet's sample names,
          # which is what we want
          x <- read.flowSet(files=match$filename, path=path)
          spillmat <- spillover_match(x, path=path, 
                                   fsc=fsc, ssc=ssc, 
                                   matchfile=matchfile)
          })



## ==========================================================================
## plot method: We actually need to attach flowViz to do the plotting
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("plot",
		signature=signature(x="flowSet",
				y="ANY"),
		definition=function(x, y, ...)
		{
			message("For plotting, please attach the 'flowViz' package.\n",
					"   i.e., 'library(flowViz)'")
		})


## ==========================================================================
## Set and replace the identifier from the environment
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("identifier",
		signature=signature(object="flowSet"),
		definition=function (object)
		{
			if(!"_.name._" %in% ls(object@frames))
				"anonymous"
			else
				object@frames[["_.name._"]]
		})

#' @export
setReplaceMethod("identifier",
		signature=signature(object="flowSet"),
		definition=function (object, value) 
		{
			object@frames[["_.name._"]] <- value
			object
		})


## ==========================================================================
## Normalize a flowSet using a normalization object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("normalize",
		signature=signature(data="flowSet",
				x="normalization"),
		definition=function (data, x)
		{
			parms <- parameters(x)
			args <-  x@arguments
			args$x <- copyFlowSet(data)
			args$parameters <- parms
			do.call(x@normFunction, args)
		})
