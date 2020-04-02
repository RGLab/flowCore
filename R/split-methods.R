## ==========================================================================
## The splitting operation in the context of 'flowFrames' and 'flowSets' is
## the logical extension of subsetting. While the latter only returns
## the events contained within a gate, the former splits the data into
## the groups of events contained within and those not contained within a
## particular gate. This concept is extremely useful in applications where
## gates describe the distinction between positivity and negativity for a
## particular marker.
## Splitting has a special meaning for gates that result in 
## 'multipleFilterResults', in which case simple subsetting doesn't make
## much sense (there are multiple populations that are defined by the gate
## and it is not clear which of those should be used for the subsetting
## operation). Accordingly, splitting of multipleFilterResults creates
## multiple subsets. The argument 'population' can be used to limit the
## output to only one or some of the resulting subsets. It takes as values
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


#' Methods to split flowFrames and flowSets according to filters
#' 
#' 
#' Divide a flow cytometry data set into several subset according to the
#' results of a filtering operation. There are also methods available to split
#' according to a factor variable.
#' 
#' 
#' The splitting operation in the context of \code{\linkS4class{flowFrame}}s
#' and \code{\linkS4class{flowSet}}s is the logical extension of subsetting.
#' While the latter only returns the events contained within a gate, the former
#' splits the data into the groups of events contained within and those not
#' contained within a particular gate. This concept is extremely useful in
#' applications where gates describe the distinction between positivity and
#' negativity for a particular marker.
#' 
#' The flow data structures in \code{flowCore} can be split into subsets on
#' various levels:
#' 
#' \code{\linkS4class{flowFrame}}: row-wise splitting of the raw data. In most
#' cases, this will be done according to the outcome of a filtering operation,
#' either using a filter that identifiers more than one sub-population or by a
#' logical filter, in which case the data is split into two populations: "in
#' the filter" and "not in the filter". In addition, the data can be split
#' according to a factor (or a numeric or character vector that can be coerced
#' into a factor).
#' 
#' \code{\linkS4class{flowSet}}: can be either split into subsets of
#' \code{\linkS4class{flowFrame}}s according to a factor or a vector that can
#' be coerced into a factor, or each individual \code{\linkS4class{flowFrame}}
#' into subpopulations based on the \code{\linkS4class{filter}}s or
#' \code{\linkS4class{filterResult}}s provided as a list of equal length.
#' 
#' Splitting has a special meaning for filters that result in
#' \code{\linkS4class{multipleFilterResult}}s or
#' \code{\linkS4class{manyFilterResult}}s, in which case simple subsetting
#' doesn't make much sense (there are multiple populations that are defined by
#' the gate and it is not clear which of those should be used for the
#' subsetting operation). Accordingly, splitting of multipleFilterResults
#' creates multiple subsets. The argument \code{population} can be used to
#' limit the output to only one or some of the resulting subsets. It takes as
#' values a character vector of names of the populations of interest. See the
#' documentation of the different filter classes on how population names can be
#' defined and the respective default values. For splitting of
#' \code{\linkS4class{logicalFilterResult}}s, the \code{population} argument
#' can be used to set the population names since there is no reasonable default
#' other than the name of the gate. The content of the argument \code{prefix}
#' will be prepended to the population names and '+' or '-' are finally
#' appended allowing for more flexible naming schemes.
#' 
#' The default return value for any of the \code{split} methods is a list, but
#' the optional logical argument \code{flowSet} can be used to return a
#' \code{\linkS4class{flowSet}} instead. This only applies when splitting
#' \code{\linkS4class{flowFrame}}s, splitting of \code{\linkS4class{flowSet}}s
#' always results in lists of \code{\linkS4class{flowSet}} objects.
#' 
#' @name split-methods
#' @aliases split-methods split split,flowFrame,ANY-method
#' split,flowFrame,factor-method split,flowFrame,character-method
#' split,flowFrame,numeric-method split,flowFrame,filter-method
#' split,flowFrame,logicalFilterResult-method
#' split,flowFrame,manyFilterResult-method
#' split,flowFrame,multipleFilterResult-method split,flowSet,ANY-method
#' split,flowSet,character-method split,flowSet,factor-method
#' split,flowSet,list-method split,flowSet,numeric-method
#' split,flowSet,filter-method split,flowSet,filterResult-method
#' @docType methods
#' @usage NULL
#' @section Methods:
#' 
#' \code{\link{flowFrame}} methods:
#' 
#' \describe{
#' 
#' \item{split(x = "flowFrame", f = "ANY", drop = "ANY")}{ Catch all input and cast an
#' error if there is no method for \code{f} to dispatch to. }
#' 
#' \item{split(x = "flowFrame", f = "factor", drop = "ANY")}{ Split a
#' \code{\link{flowFrame}} by a factor variable. Length of \code{f} should be
#' the same as \code{nrow(x)}, otherwise it will be recycled, possibly leading
#' to undesired outcomes. The optional argument \code{drop} works in the usual
#' way, in that it removes empty levels from the factor before splitting.}
#' 
#' \item{split(x = "flowFrame", f = "character", drop = "ANY")}{ Coerce \code{f} to a
#' factor and split on that. }
#' 
#' \item{split(x = "flowFrame", f = "numeric", drop = "ANY")}{ Coerce \code{f} to a
#' factor and split on that. }
#' 
#' \item{split(x = "flowFrame", f = "filter", drop = "ANY")}{ First applies the
#' \code{\linkS4class{filter}} to the \code{\linkS4class{flowFrame}} and then
#' splits on the resulting \code{\linkS4class{filterResult}} object. }
#' 
#' \item{split(x = "flowFrame", f = "logicalFilterResult", drop = "ANY")}{ Split into
#' the two subpopulations (in and out of the gate). The optional argument
#' \code{population} can be used to control the names of the results. }
#' 
#' \item{split(x = "flowFrame", f = "manyFilterResult", drop = "ANY")}{ Split into the
#' several subpopulations identified by the filtering operation. Instead of
#' returning a list, the additional logical argument codeflowSet makes the
#' method return an object of class \code{\linkS4class{flowSet}}. The optional
#' \code{population} argument takes a character vector indicating the
#' subpopulations to use for splitting (as identified by the population name in
#' the \code{filterDetails} slot).}
#' 
#' \item{split(x = "flowFrame", f = "multipleFilterResult", drop = "ANY")}{ Split into
#' the several subpopulations identified by the filtering operation. Instead of
#' returning a list, the additional logical argument codeflowSet makes the
#' method return an object of class \code{\linkS4class{flowSet}}. The optional
#' \code{population} argument takes a character vector indicating the
#' subpopulations to use for splitting (as identified by the population name in
#' the \code{filterDetails} slot). Alternatively, this can be a list of
#' characters, in which case the populations for each list item are collapsed
#' into one \code{\linkS4class{flowFrame}}.}
#' 
#' }
#' 
#' \code{\linkS4class{flowSet}} methods:
#' 
#' \describe{
#' 
#' \item{split(x = "flowSet", f = "ANY", drop = "ANY")}{ Catch all input and cast an
#' error if there is no method for \code{f} to dispatch to.  }
#' 
#' \item{split(x = "flowSet", f = "factor", drop = "ANY")}{ Split a
#' \code{\link{flowSet}} by a factor variable. Length of \code{f} needs to be
#' the same as \code{length(x)}. The optional argument \code{drop} works in the
#' usual way, in that it removes empty levels from the factor before splitting.
#' }
#' 
#' \item{split(x = "flowSet", f = "character", drop = "ANY")}{ Coerce \code{f} to a
#' factor and split on that. }
#' 
#' \item{split(x = "flowSet", f = "numeric", drop = "ANY")}{ Coerce \code{f} to a
#' factor and split on that. }
#' 
#' \item{split(x = "flowSet", f = "list", drop = "ANY")}{ Split a
#' \code{\link{flowSet}} by a list of \code{\linkS4class{filterResult}}s (as
#' typically returned by filtering operations on a
#' \code{\linkS4class{flowSet}}). The length of the list has to be equal to the
#' length of the \code{\linkS4class{flowSet}} and every list item needs to be a
#' \code{\linkS4class{filterResult}} of equal class with the same parameters.
#' Instead of returning a list, the additional logical argument codeflowSet
#' makes the method return an object of class \code{\linkS4class{flowSet}}. The
#' optional \code{population} argument takes a character vector indicating the
#' subpopulations to use for splitting (as identified by the population name in
#' the \code{filterDetails} slot).  Alternatively, this can be a list of
#' characters, in which case the populations for each list item are collapsed
#' into one \code{\linkS4class{flowFrame}}. Note that using the
#' \code{population} argument implies common population names for
#' all\code{\linkS4class{filterResult}}s in the list and there will be an error
#' if this is not the case. }
#' 
#' }
#' @author F Hahne, B. Ellis, N. Le Meur
#' @keywords methods
#' @examples
#' 
#' data(GvHD)
#' qGate <- quadGate(filterId="qg", "FSC-H"=200, "SSC-H"=400)
#' 
#' ## split a flowFrame by a filter that creates
#' ## a multipleFilterResult
#' samp <- GvHD[[1]]
#' fres <- filter(samp, qGate)
#' split(samp, qGate)
#' 
#' ## return a flowSet rather than a list
#' split(samp, fres, flowSet=TRUE)
#' 
#' ## only keep one population
#' names(fres)
#' ##split(samp, fres, population="FSC-Height+SSC-Height+")
#' 
#' 
#' ## split the whole set, only keep two populations
#' ##split(GvHD, qGate, population=c("FSC-Height+SSC-Height+",
#' ##"FSC-Height-SSC-Height+"))
#' 
#' ## now split the flowSet according to a factor
#' split(GvHD, pData(GvHD)$Patient)
#' 
#' @export




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


## Split on logicalFilterResults. This will divide the data set into those
## events that are contained within the gate and those that are not. 
#' @export
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
          keyword(out[[1]])$GUID <-
              sprintf("%s (%s)", identifier(out[[1]]), names(f)[1])
          keyword(out[[2]])$GUID <-
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
#' @export
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
              keyword(combined)$GUID <-
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
#' @export
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
#' @export
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

#' @export
setMethod("split",
          signature=signature(x="flowFrame",
                              f="numeric"),
          definition=function(x, f, drop=FALSE, prefix=NULL,
                              flowSet=FALSE, ...)
          split(x, factor(f)))

#' @export
setMethod("split",
          signature=signature(x="flowFrame",
                              f="character"),
          definition=function(x, f, drop=FALSE, prefix=NULL,
                              flowSet=FALSE, ...)
          split(x, factor(f)))


## Everything else should stop with an error
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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

#' @export
setMethod("split",
          signature=signature(x="flowSet",
                              f="numeric"),
          definition=function(x, f, drop=FALSE, ...)
          split(x, factor(f)))

#' @export
setMethod("split",
          signature=signature(x="flowSet",
                              f="character"),
          definition=function(x, f, drop=FALSE, ...)
          split(x, factor(f)))

#' @export
setMethod("split",
          signature=signature(x="flowSet",
                              f="filterResult"),
          definition=function(x, f, drop=FALSE, ...)
          stop("Can't split a flowSet by a single filterResult.\n",
               "Either provide list of filterResults or a single filter.",
               call.=FALSE))

