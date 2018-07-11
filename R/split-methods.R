## ==========================================================================
## The splitting operation in the context of 'flowFrames' and 'flowSets' is
## the logical extension of subsetting. While the latter only returns
## the events contained within a gate, the former splits the data into
## the groups of events cotained within and those not contained within a
## particular gate. This concept is extremely useful in applications where
## gates describe the distinction between positivity and negativity for a
## particual marker.
## Splitting has a special meaning for gates that result in 
## 'multipleFilterResults', in which case simple subsetting doesn't make
## much sense (there are multiple populations that are defined by the gate
## and it is not clear which of those should be used for the subsetting
## operation). Accordingly, splitting of multipleFilterResults creates
## multiple subsets. The argument 'population' can be used to limit the
## ouput to only one or some of the resulting subsets. It takes as values
## a character vector of names of the populations of interest. See the
## documentation of the different filter classes on how population names
## can be defined and the respective default values. For splitting of
## 'logicalFilterResults', the 'population' argument can be used to set
## the population names since there is not reasonable default other than the
## name of the gate. The content of the argument 'prefix' will be prepended
## to the population names and '+' or '-' are finally appended allowing for
## more flexible naming schemes.
## Further control of the output is provided by the argument 'flowSet',
## which defines whether individual subsets should be returned in the form
## of a list (the default) or whether they should be coerced into objects
## of class 'flowSet'. This only applies when splitting 'flowFrames',
## splitting of 'flowSets' always results in lists of 'flowSet' objects.
## ==========================================================================






## ==========================================================================
## split methods for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We actually split on filterResults and multipleFilterResults, so filters
## have to be evaluated first. Note that the 'drop' argument is silently
## ignored when splitting by filter since it doesn't have a clear meaning in
## this application. It does have the expected behaviour when splitting by
## factors, though.

## Evaluate the filter first and split on the filterResult
setMethod("split",
          signature=signature(x="flowFrame",
                              f="filter"),
          definition=function(x, f, drop=FALSE, ...)
          split(x, filter(x, f), drop=drop, ...))

## Evalute the filterSet first and split on the filterResult
## FIXME: Is that really what we want? And what is the final output of
## a filterSet filtering operation anyways? Just the leaves???
setMethod("split",
          signature=signature(x="flowFrame",
                              f="filterSet"),
          definition=function(x, f, drop=FALSE, ...)
          split(x, filter(x, f), drop=drop, ...))

## Split on logicalFilterResults. This will divide the data set into those
## events that are contained within the gate and those that are not. 
setMethod("split",
          signature=signature(x="flowFrame",
                              f="logicalFilterResult"),
          definition=function(x, f, drop=FALSE, population=NULL, prefix=NULL,
                              flowSet=FALSE, ...)
      {
          ## take filterID as default population name and prepend prefix
          ## if necessary
          if(is.null(population))
              population <- paste(f@filterId, c("+","-"), sep="")
          population <- unlist(population)
          if(!is.character(population))
              stop("'population' must be character scalar.", call.=FALSE)
          if(length(population)>2)
              stop("Population argument when splitting on logicalFilterResults",
                   " must be of length 2", call.=FALSE)
          allThere <- population %in% names(f)
          if(!all(allThere))
              stop("The following are not valid population names in this ",
                   "filterResult:\n\t", paste(population[!allThere],
                                        collapse="\n\t"), call.=FALSE)
          if(!is.null(prefix)){
              if(!is.character(prefix))
                  stop("'prefix' must be character vector.", call.=FALSE) 
              population <- paste(prefix[1:(min(2, length(prefix)))],
                                  population, sep="")
          }
          out <- structure(list(x[f@subSet, ], x[!f@subSet, ]),
                           names=names(f))
          description(out[[1]])$GUID <-
              sprintf("%s (%s)", identifier(out[[1]]), names(f)[1])
          description(out[[2]])$GUID <-
              sprintf("%s (%s)", identifier(out[[2]]), names(f)[2])
          keep <- match(population, names(f))
          out <- out[keep]
          if(flowSet){
              out <- flowSet(out)
              phenoData(out)$population <- population
              sampleNames(out) <- population
              varMetadata(out)["population", "labelDescription"] <-
                  "population identifier produced by splitting"
          }
          return(out)
      })

## Split on multipleFilterResults. The argument 'population' can be used to
## select only certain subpopulations
setMethod("split",
          signature=signature(x="flowFrame",
                              f="multipleFilterResult"),
          definition=function(x, f, drop=FALSE, prefix=NULL, flowSet=FALSE,
                             population=NULL, ...)
      {
          if(is.null(population))
              population <- names(f)
          else if(!all(sapply(population, is, "character")))
              stop("'population' must be a single character vector ",
                   "or a list of character vectors", call.=FALSE)
          if(!is.list(population)){
              n <- population
              population <- as.list(population)
              names(population) <- n
          }
          pop <- unique(unlist(population))
          allThere <- pop %in% names(f)
          if(!all(allThere))
              stop("The following are not valid population names in this ",
                   "filterResult:\n\t", paste(pop[!allThere],
                                        collapse="\n\t"), call.=FALSE)
          np <- length(population)
          if(is.null(prefix))
              nn <- names(population)
          else{
              if(!is.character(prefix))
                  stop("'prefix' must be character vector.", call.=FALSE)
              prefix <- rep(prefix, np)
              nn <- paste(prefix[1:np], names(population), sep="")
          }
          out <- vector(np, mode="list")
          names(out) <- nn
          i <- 1
          for(p in population){        
              tmp <- lapply(p, function(i) x[f[[i]], ])
              combined <- as(as(tmp, "flowSet"), "flowFrame")
              cn <- match("Original", colnames(combined))
              if(!is.na(cn))
                  combined <- combined[,-cn]
              description(combined)$GUID <-
                  sprintf("%s (%s)", identifier(tmp[[1]]),
                          paste(p, collapse=","))
              out[[i]] <- combined
              i <- i+1
          }
          if(flowSet){
              out <- flowSet(out)
              phenoData(out)$population <- names(population)
              sampleNames(out) <- names(population)
              varMetadata(out)["population", "labelDescription"] <-
                  "population identifier produced by splitting"
          }
          return(out)
      })


## Split on manyFilterResults. FIXME: Need to take a closer look at this
setMethod("split",
          signature=signature(x="flowFrame",
                              f="manyFilterResult"),
          definition=function(x, f, drop=FALSE, prefix=NULL,
                              flowSet=FALSE, population=NULL, ...)
      {
          m <- getMethod("split", signature=signature(x="flowFrame",
                                                      f="multipleFilterResult"))
          m(x=x, f=f, drop=drop, prefix=prefix, flowSet=flowSet,
            population=population, ...)   
      })


## Split on a factor, or on a vector that is easily coerced into a factor.
## This is just for completeness.
setMethod("split",
			 signature = signature(x = "flowFrame",
			 							 f = "factor"),
			 definition = function(x,
			 							 f,
			 							 drop = FALSE,
			 							 prefix = NULL,
			 							 flowSet = FALSE,
			 							 ...
			 )
      {      
          if (drop)
              f <- factor(f)
          nn  <- levels(f)
          out <- structure(lapply(nn,
          								function(i) x[f == i,]),
                           names = if(is.null(prefix)) nn else
                           paste(prefix, nn, sep = ""))
          if(flowSet) {
              print(data.frame(name=nn, split=seq_along(nn), row.names=nn))
              flowSet(out,phenoData=new("AnnotatedDataFrame",
                          data=data.frame(name=I(nn), split=seq_along(nn),
                          row.names=nn),
                          varMetadata=data.frame(labelDescription=I(c("Name",
                                                 "Split")),
                          row.names=c("name","split"))))
          } else out
      })

setMethod("split",
          signature=signature(x="flowFrame",
                              f="numeric"),
          definition=function(x, f, drop=FALSE, prefix=NULL,
                              flowSet=FALSE, ...)
          split(x, factor(f)))

setMethod("split",
          signature=signature(x="flowFrame",
                              f="character"),
          definition=function(x, f, drop=FALSE, prefix=NULL,
                              flowSet=FALSE, ...)
          split(x, factor(f)))


## Everything else should stop with an error
setMethod("split",
          signature=signature(x="flowFrame", f="ANY"),
          definition=function(x, f, drop=FALSE, prefix=NULL,...) 
          stop("Don't know how to split a 'flowFrame' by an object of class '",
               class(f), "'.", call.=FALSE))







## ==========================================================================
## split methods for flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## By default, we try to split the flowSet by using fsApply. As a special
## case, we need to be able to split according to a list of filtes or
## filterResults, but we also need to make sure that this list is valid, i.e.,
## only contains filters or filterResults of the same class, using the same
## parameters, etc.
## Splitting a flowSet always returns a list of flowSets.
## FIXME: How do we treat the cases in which multipleFilterResult produce
## different numbers of populations? We can't collapse to a flowSet any more,
## do we want lists of lists? Or should this be disallowed completely?


## Try and split a flowSet by whatever comes your way...
setMethod("split",
          signature=signature(x="flowSet",
                              f="ANY"),
          definition=function(x, f, drop=FALSE, population=NULL,
                             prefix=NULL, ...)
      {
          ## Split always returns a list
          sample.name <- sampleNames(x)
          fsApply(x,function(y) {
              l <- split(y, f, drop, population, prefix, ...)
              names(l) <- paste(names(l), "in", sample.name[1])
              sample.name <<- sample.name[-1]
              l
          }, simplify=FALSE)
      })


## Split a flowSet by a single filter, by first creating a list of
## filterResult and then working our way through that in the next
## method.
setMethod("split",
          signature=signature(x="flowSet",
                              f="filter"),
          definition=function(x, f, drop=FALSE, population=NULL,
                              prefix=NULL, ...)
      {
         fres <- filter(x,f)
         split(x, fres, population=population, prefix=prefix,
               ...)
      })

## Check if two filters are compatible, i.e. are they of the same class,
## do they share the same parameters...
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
            stop("Filter parameters don't match:\n",
                 paste(parameters(ff1), collapse=", "), " vs. ",
                 paste(parameters(ff2), collapse=", "), call.=FALSE)
    }
}

## Split a flowSet according to a list of filters or filterResults
## of equal length, We make sure that this list makes sense.
## FIXME: Eventually, this function should be deprecated since we represent
## filterResult for a flowSet in filterResultLists now.
setMethod("split",
          signature=signature(x="flowSet",
                              f="list"),
          function(x, f, drop=FALSE, population=NULL,
                   prefix=NULL, ...)
      {
          ## A lot of sanity checking up front
          sample.name <- sampleNames(x)
          lf <- length(f)
          lx <- length(x)
          if(lf!=lx)
              stop("list of filterResults or filters must be same ",
                   "length as flowSet.", call.=FALSE)
          if(!all(sapply(f, is, "filter")))
              stop("Second argument must be list of filterResults or filters,",
                   call.=FALSE)
          lapply(f, compatibleFilters,  f[[1]])
          ## split everything or just some populations
          ## (if multipleFilterResult)
          if(is.null(population)){
              if(!is.null(names(f[[1]])))
                  population <- names(f[[1]])
              else
                  population <- c("positive", "negative")
          } else if(!all(sapply(population, is, "character")))
              stop("'population' must be a single character vector ",
                   "or a list of character vectors", call.=FALSE)
          if(!is.list(population)){
              n <- population
              population <- as.list(population)
              names(population) <- n
          }
          ## FIXME: Do we want to allow for different names when splitting
          ## flowSets by multipleFilterResults?
          if(lf>1 && !identical(unique(as.vector(sapply(f, names))),
                                names(f[[1]]))){
              for(i in 2:lf)
                  names(f[[i]]) <- names(f[[1]])
              warning("Filtering operation produced non-unique population ",
                      "names.\n  Using names of the first frame now.\n",
                      "  Please check parameter descriptions in the ",
                      "parameter slots\n  of the individual flowFrames.",
                      call.=FALSE)
          }
          finalRes <- vector(mode="list", length=length(population))
          names(finalRes) <- names(population)
          for(p in seq_along(population)){
              tp <- population[p]
              res <- vector(mode="list", length=lf)
              for(i in 1:lf){
                  l <- split(x[[i]], f[[i]], population=tp,
                             prefix=prefix, flowSet=FALSE, ...)
                  res[[i]] <- l[[1]]
                  if(!is.null(prefix)){
                      if(is.logical(prefix) && prefix)
                          names(res)[i] <- paste(names(l), "in", sample.name[i])
                      else if(is.character(prefix))
                          names(res)[i] <- paste(prefix, sample.name[i])
                  }else
                  names(res)[i] <- sample.name[i]   
              }
              np <- names(population)[p]
              finalRes[[np]] <- flowSet(res, phenoData=phenoData(x))
              phenoData(finalRes[[np]])$population <- np
              varMetadata(finalRes[[np]])["population", "labelDescription"] <-
                  "population identifier produced by splitting"
          }
          return(finalRes)
      })

## FIXME: This should replace the above list method completely at some point
##  NOTE: This method now replaces the one with signature("flowSet","list") [Greg Finak <gfinak@fhcrc.org>, 04/13/2011]
setMethod("split",
          signature=signature(x="flowSet",
                              f="filterResultList"),
          definition=function(x, f, drop=FALSE, population=NULL,
                              prefix=NULL, ...)
      {
		lf<-length(f)
		sample.name <- sampleNames(x)
		if(length(x)!=length(f)){
			stop("filterResultList and flowSet must be same ",
                   "length.", call.=FALSE)		
		}
		lapply(f, flowCore:::compatibleFilters,  f[[1]])
          if(is.null(population)){
              if(all(unlist(lapply(f,function(q)!is.null(names(q))))))
 				population<-unique(unlist(lapply(f,names)))
              else
                  population <- c("positive", "negative")
          } else if(!all(sapply(population, is, "character")))
              stop("'population' must be a single character vector ",
                   "or a list of character vectors", call.=FALSE)
          if(!is.list(population)){
              n <- population
              population <- as.list(population)
              names(population) <- n
          }
          finalRes <- vector(mode="list", length=length(population))
          names(finalRes) <- names(population)
          for(p in seq_along(population)){
	          tp <- population[p]
              res <- vector(mode="list", length=lf)
    		for(i in 1:lf){
				l <- try(split(x[[i]], f[[i]], population=tp,
                           prefix=prefix, flowSet=FALSE, ...),silent=TRUE)
				if(inherits(l,"try-error")){
					if(geterrmessage()==paste("Error : The following are not valid population names in this filterResult:\n\t",tp,"\n",sep="")){
						message("Creating an empty flowFrame for population ",tp,"\n")
						#Create an empty flowFrame
						l<-x[[i]][0,];
						identifier(l)<-paste(identifier(l),paste("(",tp,")",sep=""),sep=" ")
						l<-list(l);
					}else
						stop("Can't split flowFrame ",sampleNames(x[i])," on population ",tp);
				}
                res[[i]] <- l[[1]]
                if(!is.null(prefix)){
                    if(is.logical(prefix) && prefix)
                        names(res)[i] <- paste(names(l), "in", sample.name[i])
                    else if(is.character(prefix))
                        names(res)[i] <- paste(prefix, sample.name[i])
                }else
                names(res)[i] <- sample.name[i]      
			}
			np <- names(population)[p]
            finalRes[[np]] <- flowSet(res, phenoData=phenoData(x))
            phenoData(finalRes[[np]])$population <- np
            varMetadata(finalRes[[np]])["population", "labelDescription"] <-
                "population identifier produced by splitting"
		  }
          #n <- f@frameId
          #f <- f@.Data
          #names(f) <- n
          #split(x, f, drop=drop, population=NULL, prefix=NULL, ...)
		return(finalRes);
      })

## Split by frames of flowSet according to a factor, character or numeric.
## Those have to be of the same length as the flowSet. We can't allow for
## drop=TRUE, because this would create invalid sets.
setMethod("split",
          signature=signature(x="flowSet",
                              f="factor"),
          definition=function(x, f, drop=FALSE, ...)
      {
          if(!is.atomic(f) || length(f)!=length(x))
              stop("split factor must be same length as flowSet",
                   call.=FALSE) 
          gind <- split(1:length(f), f, drop=TRUE)
          res <- vector(mode="list", length=length(gind))
          for(g in seq_along(gind)){
              res[[g]] <- x[gind[[g]]]
              phenoData(res[[g]])$split <- levels(f)[g]
              varMetadata(res[[g]])["split", "labelDescription"] <-
                  "Split"
          }
          names(res) <- names(gind)
          return(res)
      })

setMethod("split",
          signature=signature(x="flowSet",
                              f="numeric"),
          definition=function(x, f, drop=FALSE, ...)
          split(x, factor(f)))

setMethod("split",
          signature=signature(x="flowSet",
                              f="character"),
          definition=function(x, f, drop=FALSE, ...)
          split(x, factor(f)))

setMethod("split",
          signature=signature(x="flowSet",
                              f="filterResult"),
          definition=function(x, f, drop=FALSE, ...)
          stop("Can't split a flowSet by a single filterResult.\n",
               "Either provide list of filterResults or a single filter.",
               call.=FALSE))

