## ==========================================================================
## split methods for flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We actually split on filterResults and multipleFilterResults
setMethod("split",
          signature("flowFrame", "filter"),
          function(x, f, drop=FALSE, ...)
          split(x, filter(x, f), drop, ...))

setMethod("split",
          signature("flowFrame","filterSet"),
          function(x, f, drop=FALSE, ...)
          split(x, filter(x, f), drop, ...))

## split on logicalFilterResults
setMethod("split",
          signature("flowFrame", "logicalFilterResult"),
          function(x, f, drop=FALSE, population=NULL, prefix=NULL,
                   flowSet=FALSE, ...)
      {
          if(is.null(population))
              population <- f@filterId
          if(!is.null(prefix))
              population <- paste(prefix,population, sep="")
          out <- structure(list(x[f@subSet, ], x[!f@subSet, ]),
                           names=c(paste(population, "+", sep=""),
                           paste(population,"-", sep="")))
          if(length(flowSet) > 0 && flowSet)
              out <- flowSet(out)
          return(if(is.list(out) && length(out)==1) out[[1]] else out)
      })

## split on multipleFilterResults, argument population can be used to
## select only certain subpopulations
setMethod("split",
          signature("flowFrame", "multipleFilterResult"),
          function(x, f, drop=FALSE, prefix=NULL, flowSet=FALSE,
                   population=NULL, ...)
      {
          if(is.null(population))
              population <- names(f)
          else if(!all(population %in% names(f)))
              stop("Population(s) not valid in this filter", call.=FALSE)
          if(is.null(prefix))
              nn <- population
          else
              nn <- paste(prefix, population, sep="")
          tmp <- lapply(population, function(i) x[f[[i]], ])
          out <- structure(tmp, names=nn)
          if(length(flowSet) > 0 && flowSet)
              out <- flowSet(out)
          return(if(is.list(out) && length(out)==1) out[[1]] else out)
      })


## Split on manyFilterResults. FIXME: Need to take a closer look at this
setMethod("split", signature("flowFrame","manyFilterResult"),
          function(x, f, drop=FALSE, prefix=NULL, flowSet=FALSE, ...)
      {
          ##If drop is TRUE then only use results without children
          if(drop)
              nn <- rownames(f@dependency)[rowSums(f@dependency)==0]
          else
              nn <- names(f)
          out <- structure(lapply(nn, function(i) x[f[[i]], ]),
                           names=if(is.null(prefix)) nn else
                           paste(prefix,nn,sep=""))
          if(length(flowSet) > 0 && flowSet) {
              print(data.frame(name=nn, as.data.frame(f)[nn, ], row.names=nn))
              out <- flowSet(out, phenoData=new("AnnotatedDataFrame",
                                  data=data.frame(name=I(nn),
                                  as.data.frame(f)[nn, ],
                                  row.names=nn),
                                  varMetadata=data.frame(labelDescription=I(
                                                         c("Name", "Filter")),
                                  row.names=c("name", "filter"))))
          }
          return(if(is.list(out) && length(out)==1) out[[1]] else out)
      })


## Split on a factor, or on a vector that is easily coerced into a factor.
## This is just for completeness.
setMethod("split", signature("flowFrame", "factor"),
          function(x, f, drop=FALSE, prefix=NULL, flowSet=FALSE, ...)
      {      
          nn  <- levels(f)
          out <- structure(lapply(nn,function(i) x[f==i,]),
                           names=if(is.null(prefix)) nn else
                           paste(prefix, i, sep=""))
          if(length(flowSet) > 0 && flowSet) {
              print(data.frame(name=nn, split=seq_along(nn), row.names=nn))
              flowSet(out,phenoData=new("AnnotatedDataFrame",
                          data=data.frame(name=I(nn), split=seq_along(nn),
                          row.names=nn),
                          varMetadata=data.frame(labelDescription=I(c("Name",
                                                 "Split")),
                          row.names=c("name","split"))))
          } else out
      })

setMethod("split",signature("flowFrame","numeric"),
          function(x, f, drop=FALSE, prefix=NULL, flowSet=FALSE, ...)
          split(x, factor(f)))

setMethod("split",signature("flowFrame","character"),
          function(x, f, drop=FALSE, prefix=NULL, flowSet=FALSE, ...)
          split(x, factor(f)))


## Everything else should stop with an error
setMethod("split",
          signature("flowFrame","ANY"),
          function(x, f, drop=FALSE, prefix=NULL,...) 
          stop("invalid type for flowFrame split"))










## ==========================================================================
## split methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## split a flowSet by whatever comes your way...
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


## split a flowSet by a single filter, by first creating a filterResult
setMethod("split",signature("flowSet","filter"),
          function(x, f, drop=FALSE, population=NULL, prefix=NULL,
                   flowSet=FALSE, ...)
      {
         fres <- filter(x,f)
         split(x, fres, population=population, prefix=prefix, flowSet=flowSet,
               ...)
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
                             prefix=prefix, flowSet=FALSE, ...)
                  res[[i]] <- l
                  if(!is.null(prefix)){
                      if(is.logical(prefix) && prefix)
                          names(res)[i] <- paste(names(l), "in", sample.name[i])
                      else if(is.character(prefix))
                          names(res)[i] <- paste(prefix, sample.name[i])
                  }else
                  names(res)[i] <- sample.name[i]   
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
          gind <- split(1:length(f), f, drop=TRUE)
          res <- vector(mode="list", length=length(gind))
          for(g in seq_along(gind))
              res[[g]] <- x[gind[[g]]]
          names(res) <- names(gind)
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
