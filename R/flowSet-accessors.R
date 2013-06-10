## ==========================================================================
## flowSets are basically lists flowFrames
## ==========================================================================






## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## to flowSet
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
setMethod("colnames",
		signature=signature(x="flowSet"),
		definition=function(x, do.NULL="missing", prefix="missing")
			x@colnames)

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



## ==========================================================================
## Allow for the extraction and replacement of phenoData
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("phenoData",
		signature=signature(object="flowSet"),
		definition=function(object) object@phenoData)

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
			object@phenoData <- value
			object
		})



## ==========================================================================
## directly access the pData data frame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("pData",
		signature=signature(object="flowSet"),
		definition=function(object) pData(object@phenoData))

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
setMethod("varLabels",
		signature=signature(object="flowSet"),
		function(object) varLabels(phenoData(object)))

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

setMethod("varMetadata",
		signature=signature(object="flowSet"),
		definition=function(object) varMetadata(phenoData(object)))

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
setMethod("sampleNames",
		signature=signature(object="flowSet"),
		definition=function(object) 
			sampleNames(phenoData(object)))

## Note that the replacement method also replaces the GUID for each flowFrame
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
                        pd$name <- value
			object@phenoData <- pd
			object@frames <- env
			return(object)
		})


## ==========================================================================
## keyword method for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
								y <- as(x[[n]],"flowFrame")
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
setMethod("compensate",
		signature=signature(x="flowSet",
				spillover="ANY"),
		definition=function(x, spillover)
			fsApply(x, compensate, spillover))



## ==========================================================================
## Transformation methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
		signature=signature(`_data`="flowSet"),
		definition=function(`_data`,...)
		{
			fsApply(`_data`,transform,...)
		})

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
## These methods apply single filters, filterSets or lists of filters to
## flowSet object. In all cases, the output of the filtering operation is
## a filterResultList
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## for filters
setMethod("filter",
		signature=signature(x="flowSet",
				filter="filter"),
		definition=function(x, filter)
		{
			if(!all(parameters(filter) %in% colnames(x)))
				stop("parameters in the filter definition don't ",
						"match the parameters in the flowSet", call.=FALSE)
			res <- fsApply(x,function(x) filter(x,filter))
			return(new("filterResultList", .Data=res, frameId=sampleNames(x),
							filterId=identifier(filter)))
		})

## for filterSets
## FIXME: Need to check that everything still works after introduction
## of filterResultLists
setMethod("filter",
		signature=signature(x="flowSet",
				filter="filterSet"),
		definition=function(x, filter)
		{
			res <- fsApply(x, function(x) filter(x, filter))
			return(new("filterResultList", .Data=res, frameId=sampleNames(x),
							filterId=identifier(filter)))
		})


setMethod("filter",
		signature=signature(x="flowSet",
				filter="list"),
		definition=function(x, filter)
      {
          filter(x, filterList(filter))
      })

## for named lists of filters. Names of the list items have to correspond
## to sampleNames in the set. Filters in the filter list that can't be
## matched are ignored, for those that are missing, an "empty" dummy
## filterResult is produced
setMethod("filter",
		signature=signature(x="flowSet",
				filter="filterList"),
		definition=function(x, filter)
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
setMethod("rbind2",
		signature=signature(x="flowSet",
				y="missing"),
		definition=function(x, y) x)

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

setMethod("rbind2",
		signature=signature(x="flowSet",
				y="flowFrame"),
		definition=function(x,y)
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

setMethod("rbind2",
		signature=signature(x="flowFrame",
				y="flowSet"),
		definition=function(x,y) rbind2(y,x))






## ===========================================================================
## spillover method
## ---------------------------------------------------------------------------
setMethod("spillover",
          signature = signature(x = "flowSet"),
          definition = function(x, unstained = NULL, patt = NULL, fsc = "FSC-A",
                                ssc = "SSC-A", method = "median",
                                stain_match = c("intensity", "ordered", "regexpr"),
                                useNormFilt = FALSE, pregate = FALSE,
                                plot = FALSE, ...) {

            stain_match <- match.arg(stain_match)
            
            if (is.null(unstained)) {
              stop("Sorry, we don't yet support unstained cells blended ",
                   "with stained cells", call. = FALSE)
            } else {
              ## We often only want spillover for a subset of the columns
              allcols <- colnames(x)
              cols <- if (is.null(patt)) {
                allcols
              } else {
                grep(patt, allcols, value = TRUE)
              }

              ## Ignore these guys if they somehow got into cols.
              cols <- cols[!(cols %in% c(fsc, ssc))]
              
              ## There has got to be a better way of doing this...
              if (!is.numeric(unstained)) {
                unstained <- match(unstained, sampleNames(x))
                if (is.na(unstained)) {
                  stop("Baseline not in this set.", call. = FALSE)
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
              # channel names with the filenames of the compensation controls.
              if (stain_match == "intensity") {
                channel_order <- NA
              } else if (stain_match == "ordered") {
                channel_order <- seq_along(sampleNames(x))[-unstained]
              } else if (stain_match == "regexpr") {
                channel_order <- sapply(cols, grep, x = sampleNames(x), fixed = TRUE)
                if (!all(sapply(channel_order, length) == 1)) {
                  stop("Multiple stains match to a common compensation-control filename",
                       call. = FALSE)
                }
              }

              if (pregate) {
                if (any(is.na(channel_order))) {
                  stop("Cannot apply pregate without knowing ordering of channels. ",
                       "Match the channels to controls with 'ordered' or 'regexpr'.",
                       call. = FALSE)
                }
                require('flowStats')

                if (plot) {
                  oask <- devAskNewPage(TRUE)
                  on.exit(devAskNewPage(oask))
                }

                x_gated <- lapply(sort(channel_order), function(channel_i) {
                  flow_frame <- x[[channel_i]]
                  channel_name <- cols[which(channel_order == channel_i)]

                  # Applies flowStats:::rangeGate to select positive population
                  gate_filter <- rangeGate(flow_frame, stain = channel_name,
                                           inBetween = TRUE, borderQuant = 0,
                                           absolute = FALSE, ...)
                  if (plot) {
                    # Plots a kernel density for the current channel
                    plot(density(exprs(flow_frame)[, channel_name]),
                         xlab = channel_name, ylab = "Density",
                         main = paste("Compensation Control:", sampleNames(x)[channel_i]))

                    # Adds a vertical line to show gate
                    cutpoint <- c(gate_filter@min, gate_filter@max)
                    cutpoint <- cutpoint[is.finite(cutpoint)]
                    abline(v = cutpoint, col = "black", lwd = 3, lty = 2)
                  }
                  Subset(flow_frame, gate_filter)
                })
                x_gated <- x_gated[channel_order]
                names(x_gated) <- sampleNames(x)[channel_order]
                x <- rbind2(flowSet(x_gated), x[unstained])
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
              }
              rownames(inten) <- colnames(inten)
              inten
          }
      })



## ==========================================================================
## plot method: We actually need to attach flowViz to do the plotting
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
setMethod("identifier",
		signature=signature(object="flowSet"),
		definition=function (object)
		{
			if(!"_.name._" %in% ls(object@frames))
				"anonymous"
			else
				object@frames[["_.name._"]]
		})

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
