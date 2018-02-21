## ==========================================================================
## We make heavy use of coercion methods to unify the resolts that
## produced by the many different filter and transformation classes.
## ==========================================================================






## ==========================================================================
## We can convert a factor, logical or a numeric into a filterResult by
## selecting a specific filterResult type. This is done through the standard
## R coercion techniques. We use these methods to create the appropriate
## filterResults from the results of the various %in% methods.
## --------------------------------------------------------------------------
setAs(from="factor", to="filterResult", def=function(from)
      new("multipleFilterResult", filterId="", subSet=from))

setAs(from="logical", to="filterResult", def=function(from)
      new("logicalFilterResult", filterId="", subSet=from))

setAs(from="numeric", to="filterResult", def=function(from)
      new("randomFilterResult", filterId="", subSet=from))

setAs(from="matrix", to="filterResult", def=function(from)
      new("manyFilterResult", filterId="", subSet=from))





setAs(from="gateView", to="filterResult", def=function(from)
  {
      fres <- as(from@indices, "filterResult")
      ofres <- get(from@filterResult)
      filter <- filterDetails(get(from@filterResult), 1)[["filter"]]
      identifier(filter) <- from@frEntry
      filterDetails(fres, from@frEntry) <- list(filter=filter)
      identifier(fres) <- from@frEntry
      return(fres)
  })

setAs(from="list", to="filterResultList",
      def=function(from)
  {
      if(is.null(names(from)))
          stop("Can only coerce a named list to a filterResultList",
               call.=FALSE)
      new("filterResultList", .Data=from, frameId=names(from),
          filterId="default")
  })

setAs(from="filterResultList", to="list",
      def=function(from)
  {
      names(from@.Data) <- names(from)
      from@.Data
  })



## ==========================================================================
## We can also convert some filterResult types directly to logical types,
## though in general it is not possible. We provide the means for logical
## and random filter types. For the rest we cast useful error messages.
## --------------------------------------------------------------------------
setAs(from="filterResult", to="logical", def=function(from)
      stop("Unable to convert to a logical vector"))

setAs(from="logicalFilterResult", to="logical", def=function(from)
      from@subSet)

## This only makes sense under the assumption that the values in subSet
## are uniformly distributed
setAs(from="randomFilterResult", to="logical", def=function(from)
      runif(length(from@subSet)) < from@subSet)



## ==========================================================================
## Allow the coercion of resolvable filters (i.e. those derived from
## filterResult) to be composed and then converted into a logical vector.
## This allows for a lot of processing to be done simply using the filter
## results.
## --------------------------------------------------------------------------
setAs(from="filter", to="logical", def=function(from)
      stop("Only resolved filters can be converted to a logical vector."))

setAs(from="subsetFilter", to="logical", def=function(from)
      as(from@filters[[1]], "logical") & as(from@filters[[2]], "logical"))

setAs(from="intersectFilter", to="logical", def=function(from)
      apply(sapply(from@filters, as, Class="logical"), 1, all))

setAs(from="unionFilter", to="logical", def=function(from)
      apply(sapply(from@filters, as, Class="logical"), 1, any))

setAs(from="complementFilter", to="logical", def=function(from)
      !as(from@filters[[1]], "logical"))



## ==========================================================================
## Allows for the resolution of filterReferences and formulas
## --------------------------------------------------------------------------
setAs(from="filterReference", to="concreteFilter", def=function(from)
  {
      x <- from@env[[from@name]]
      if(is.null(x)) stop(paste("Unable to resolve filter reference:",
                                from@name))
      x
  })



setAs(from="formula", to="filter", def=function(from)
  {
      f <- as(from[[length(from)]], "filter")
      if(length(from) == 3 && from[[2]] != ".")
          f@filterId = as.character(from[[2]])
      f
  })

setAs(from="character", to="filter", def=function(from)
    filterReference(as(find(from, mode='S4'), "environment"), from))

setAs(from="name", to="filter", def=function(from)
      as(as.character(from), "filter"))

setAs(from="call","filter", def=function(from)
  {
      filters <- lapply(from[-1], as, Class="filter")
      eval(as.call(c(from[[1]], filters)))
  })



## ==========================================================================
## These exist primarily to support making copies of filterSet objects
## --------------------------------------------------------------------------
setAs(from="filterReference", to="call", def=function(from)
      as.symbol(from@name))

setAs(from="filter", to="call", def=function(from)
  {
      nam <- names(getSlots(class(from)))
      vals <- structure(lapply(nam, function(n) {
          v <- slot(object=from,n)
          if(is.call(v)) as.call(c(as.symbol("quote"), v)) else v
      }), names=nam)
      as.call(c(as.symbol("new"), class(from), vals))
  })

setAs(from="subsetFilter", to="call", def=function(from) {
    eval(as.call(c(as.symbol('~'), as.symbol(from@filterId),
                   as.call(c(as.symbol("%subset%"),
                             lapply(from@filters, as, Class="call"))))))
})

## Helper function for converting lists into binary call trees
binaryHelper <- function(op,l)
{
    x <- l[[1]]
    op <- as.symbol(op)
    for(i in 2:length(l))
        x <- as.call(c(op,x,l[[i]]))
    x
}

setAs(from="intersectFilter", to="call", def=function(from) {
    eval(as.call(c(as.symbol('~'), as.symbol(from@filterId),
                   binaryHelper('&',lapply(from@filters, as, Class="call")))))
})

setAs(from="unionFilter", to="call",def=function(from) {
    eval(as.call(c(as.symbol('~'), as.symbol(from@filterId),
                   binaryHelper('|', lapply(from@filters, as, Class="call")))))
})

setAs(from="complementFilter", to="call", def=function(from) {
    if(length(from@filters) > 1)
        stop("Whoops. Complements only work on one filter right now.")
    eval(as.call(c(as.symbol('~'), as.symbol(from@filterId),
                   as.call(c(as.symbol("!"), as(from@filters[[1]], "call"))))))
})



## ==========================================================================
## Lists to filterSets and vice versa
## --------------------------------------------------------------------------
setAs(from="filterSet", to="list", def=function(from)
  {
      nam <- ls(envir = from@env)
      out <- lapply(nam,function(n) as(from[[n]], "call"))
      names(out) <- nam
      out
  })

setAs(from="list", to="filterSet", def=function(from)
  {
      fs <- filterSet()
      e <- new.env()
      n <- names(from)
      if(is.null(n) || length(n)==0) n = rep("", length(from))
      for(i in 1:length(from)) {
          filter <- from[[i]]
          name <-if(!is.null(n[i]) && nzchar(n[i])) n[i] else NULL
          fs[[name]] <- eval(filter,e)
      }
      fs
  })



## ==========================================================================
## Convert an environment to a flowSet.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="environment", to="flowSet", def=function(from)
  {
      frameList <- ls(envir = from)
      isFrame <- sapply(frameList, function(f) is(get(f, envir = from),
                                                  "flowFrame"))
      if(!all(isFrame))
          warning("Some symbols are not flowFrames.",
                  "They will be ignored but left intact.")
      ## If specified, remove extraneous symbols from the environment
      ## before continuing
      frameList <- frameList[isFrame]
      new("flowSet", frames=from, colnames=colnames(from[[frameList[[1]]]]),
          phenoData=new("AnnotatedDataFrame",
                        data=data.frame(name=I(frameList), row.names=frameList),
                        varMetadata=data.frame(labelDescription="Name",
                                               row.names="name")))
  })



## ==========================================================================
## Convert a list to a flowSet by creating an environment and coerce THAT,
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="list", to="flowSet", def=function(from)
  {
      if(is.null(names(from)))
          names(from) <- paste("V", seq(1, length(from)), sep="")
        
      orig.sampleNames <- names(from)
      #this is the trick that forces the current sample names in alphabetic order by prepending serial numbers so that the original order of samples is preserved instead of being shuffled by list2env call (and within setAs(from="environment", to="flowSet") method 'ls` call reorder the list by alphabet order)
      names(from) <- paste(sprintf("%0.6d", seq_along(from)), names(from), sep="_")
      res <- as(list2env(from, new.env(hash=T, parent=emptyenv())), "flowSet")
      #by reassigning the original sample names, we are also using its side effect to overwrite the GUID keyword in flow data, which will prevent read.flowSet from
      #renaming the flowSet with GUID (which could be a design bug in itself).
      sampleNames(res) <- orig.sampleNames
      #restore name column in pData as well (since it is no longer taken care of by sampleNames<- method)
      pData(res)[["name"]] <- I(orig.sampleNames)
      res
  })



## ==========================================================================
## Convert a flowSet to a list
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="flowSet", to="list", def=function(from) {
    ret <- list()
    for(i in sampleNames(from))
        ret[[i]] <- from[[i]]
    return(ret)
})



## ==========================================================================
## Coerce a flowFrame to a flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="flowFrame", to="flowSet", def=function(from)
    flowSet(from))



## ==========================================================================
## Coerce a flowSet to a flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="flowSet", to="flowFrame", def=function(from)
  {
      if(length(from) == 1)
          from[[1]]
      else {
          ## The parameters need to match all frames
          params <- parameters(from[[1]])
          allParams <- fsApply(from, function(x)
                               as.character(pData(parameters(x))$name))
          if(!all(apply(allParams, 2, function(x) length(unique(x))==1)))
              stop("parameters must be the same for all frames")
          ## making sure we are not doing too many copies of the data
          lens <- fsApply(from, nrow)
          exp <- matrix(ncol=nrow(params)+1, nrow=sum(lens))
          colnames(exp) <- c(colnames(from), "Original")
          offset <- 1
          for(i in 1:length(from)){
              if(lens[[i]]>0){
                  rows <- offset:(offset+lens[i,]-1)
                  exp[rows, 1:nrow(params)] <- exprs(from[[i]])
                  exp[rows,"Original"] <- rep(i, lens[i,])
                  offset <- offset+lens[i,]
              }
          }
          repl <-  data.frame(name="Original", range=NA, minRange=1,
                              maxRange=length(from), stringsAsFactors=FALSE)
          rownames(repl) <- "Original"
          common <- intersect(colnames(repl), colnames(pData(params)))
          pData(params)["Original",common] <- repl[,common]
          pData(params)[,"desc"] <-
            c(as.character(pData(parameters(from[[1]]))[,"desc"]),
                                 "Original Frame")
          desc  <- list(description="Synthetic Frame",
                        sampleNames=sampleNames(from))
          new("flowFrame",exprs=exp,parameters=params,description=desc)
      }
  })



## ==========================================================================
## Coerce a filterSummary to a data.frame. This gets used by the toTable
## methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="filterSummary", to="data.frame", def=function(from)
      data.frame("true"=from@true, "false"=from@count-from@true,
                 "count"=from@count,"p"=from@p,
                 "q"=1-from@q,row.names=from@name))




## ==========================================================================
## coerce from a list to a transformList
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="list", to="transformList", def=function(from)
      new("transformList", transforms=from))




## ==========================================================================
## coerce from a transform object to characters (if possible). This needs
## to be recursive because parameters of transforms can again be transforms.
## We return NULL if we don
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## We know how to coerce a unitytransform
setAs(from="unitytransform", to="character", def=function(from){
  tmp <- from@parameters
  tmp
})

## We can't coerce the ratio transformation in any case, so we return a
## nullParameter
setAs(from="ratio", to="character", def=function(from){
  from <- new("nullParameter")
  selectMethod("coerce", c("transform", "character"))(from)
})

## Coercing a nullParameter gives us NA
setAs(from="nullParameter", to="character", def=function(from){
  tmp <- NA
  tmp
})

## recursively coerce the parameters slot
setAs(from="transform", to="character", def=function(from)
      {
        p <- parameters(from)
        if(is.character(p)){
          if(length(p)==1){
            return(p)
          }else
          return(new("nullParameter"))
        }else{
          return(sapply(p, as, "character"))
        }
      })



setAs(from="parameters", to="character", def=function(from){
      tmp <- sapply(from, as, "character")
      tmp
    })




## ==========================================================================
## coerce between gate representations
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setAs(from="ellipsoidGate", to="polygonGate", def=function(from)
  {
      parms <- parameters(from)
      ## get the ellipse lines
      center <- from@mean[parms]
      if(is.null(rownames(from@cov)))
          rownames(from@cov) <- colnames(from@cov)
      cov <- from@cov[parms, parms]
      radius <- from@distance
      chol.cov <- t(chol(cov))
      t <- seq(0, 2 * base::pi, length = 50)
      ans <- center +
          (chol.cov %*% rbind(x = radius * cos(t),
                              y = radius * sin(t)))
      ans <- as.data.frame(t(ans))
      names(ans) <- parms
      ## create a polygonGate
      g <- polygonGate(.gate=ans, filterId=identifier(from))
      #need do this to preserve the transform info of the original gate parameters 
      g@parameters <- from@parameters
      g
  })



setAs(from="rectangleGate", to="polygonGate", def=function(from)
  {
      bound <- rbind(from@min, c(from@max[1], from@min[2]), from@max,
                     c(from@min[1], from@max[2]))
      g <- polygonGate(.gate=bound, filterId=identifier(from))
      #need do this to preserve the transform info of the original gate parameters 
      g@parameters <- from@parameters
      g
  })

