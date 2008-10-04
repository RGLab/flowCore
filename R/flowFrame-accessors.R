## ==========================================================================
## flowFrames are the basic data structures for flow cytometry data
## ==========================================================================






## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## by indices or logical vectors
setMethod("[",
          signature=signature(x="flowFrame"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          if(drop)
              warning("Argument 'drop' ignored for subsetting of flowFrame")
          msg <- "Subset out of bounds"
          if(!missing(j)){
              if(is.logical(j))
                  if(max(which(j)) > ncol(x))
                      stop(msg, call.=FALSE)
              if(is.numeric(j))
                  if(max(abs(j)) > ncol(x))
                      stop(msg, call.=FALSE)
              if(is.character(j))
                  if(!all(j %in% colnames(x)))
                      stop(msg, call.=FALSE)
          }
          if(!missing(i))
               if(max(abs(i)) > nrow(x))
                   stop(msg, call.=FALSE)
          
          exprs(x) <-  switch(1+missing(i)+2*missing(j),
                          { exprs(x)[i, j, ..., drop=FALSE] },
                          { exprs(x)[ , j, ..., drop=FALSE] },
                          { exprs(x)[i,  , ..., drop=FALSE] },
                          { exprs(x)[ ,  , ..., drop=FALSE] } )
          x
      })

## by results of a filtering operation
setMethod("[",
          signature=signature(x="flowFrame",
                              i="filterResult"),
          definition=function(x,i,j,...,drop=FALSE)
      {
          if(missing(j))
              x[x %in% i,,...,drop=FALSE]
          else
              x[x %in% i,j,...,drop=FALSE]
      })

## by filter (which computes filter result first and applies that)
setMethod("[",
          signature=signature(x="flowFrame",
                              i="filter"),
          definition=function(x,i,j,...,drop=FALSE)
      {
          result <- filter(x,i)
          if(missing(j))
              x[result,,...,drop=FALSE]
          else
              x[result,j,...,drop=FALSE]
      })




## ==========================================================================
## the $-operator (subsetting with '$')
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"$.flowFrame" <- function(x, val){
    if(!val %in% colnames(x))
        stop("'", val, "' is not a valid parameter in this flowFrame\n",
             "Subscript out of bounds")
    x[,val]
}




## ==========================================================================
## accessor and replace methods for slot exprs
## Note that only numeric matrices are accepted as replacement values
## and that the method checks if the colnames match the parameter
## definition in the parameters slot. Implicit subsetting is allowed
## (i.e. less columns in the replacement value compared to the original
## flowFrame, but all have to be defined there) 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("exprs",
          signature=signature(object="flowFrame"),
          definition=function(object)
          object@exprs
          )

setReplaceMethod("exprs",
                 signature=signature(object="flowFrame",
                                     value="matrix"),
                 definition=function(object, value)
             {
                 if(!is.numeric(value))
                     stop("replacement value for the 'exprs' slot ",
                          "of 'flowFrame' objects must be numeric",
                          "matrix", call.=FALSE)
                 mt <- match(colnames(value), parameters(object)$name)
                 if(length(mt)==0)
                     stop("colnames missing in replacement value. ",
                          "Unable to match data columns", call.=FALSE)
                 if(any(is.na(mt)))
                     stop("the following parameters are not defined in the ",
                          "'flowFame' object:\n\t",
                          paste(colnames(value)[is.na(mt)], collapse=", "),
                          "\nunable to replace", call.=FALSE)
                 object@parameters <- object@parameters[mt,,drop=FALSE]
                 object@exprs <- value
                 return(object)
             })

## throw meaningful error when trying to replace with anything other than matrix
setReplaceMethod("exprs",
                 signature=signature(object="flowFrame", value="ANY"),
                 definition=function(object, value)
                 stop("Replacement value for the 'exprs' slot of a ",
                      "'flowFrame' object must be \n  a numeric matrix with ",
                      "colnames matching at least a subset of the orginal",
                      "\n  columns."))



## ==========================================================================
## accessor and replace methods for slot description
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("description",
          signature=signature(object="flowFrame"),
          function(object, hideInternal=FALSE){
              if(!hideInternal)
                  object@description
              else{
                  sel <- grep("^\\$", names(object@description))
                  object@description[-sel]
              }})

## replace description entries
descError <- "Replacement value must be a named list."
setReplaceMethod("description",
                 signature=signature(object="flowFrame",
                                     value="list"),
                 definition=function(object, value)
             {
                 n <- names(value)
                 if(length(n) == 0)
                     stop(descError, call.=FALSE)
                 object@description[n] <- value
                 return(object) })

setReplaceMethod("description",
                 signature=signature(object="flowFrame",
                                     value="ANY"),
                 definition=function(object, value)
                 stop(descError, call.=FALSE))


## ==========================================================================
## accessor methods for individual items in the description slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## select keywords by name
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="character"),
          function(object, keyword)
              structure(object@description[keyword], names=keyword)
          )

## select keywords by function
## FIXME: What is the idea behind this? Don't see the use case...
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="function"),
          definition=function(object,keyword)
          keyword(object,object@description)
          )

## select keywords by combination of name and function
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="list"),
          definition=function(object,keyword)
      {
          sapply(keyword,function(k) {
              if(is.character(k))
                  object@description[k]
              else if(is.function(k))
                  k(object,object@description,k)
              else NA
          })
      })

## this is equivalent to the description method
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="missing"),
          function(object)
              object@description
          )

## replace keywords
kwdError <- "Replacement value must be a named character vector or list."
setReplaceMethod("keyword",
                 signature=signature(object="flowFrame", value="character"),
                 definition=function(object, value)
             {
                 keyword(object) <- as.list(value)
                 return(object)
             })

setReplaceMethod("keyword",
                 signature=signature(object="flowFrame",
                                     value="list"),
                 definition=function(object, value)
             {
                 n <- names(value)
                 if(length(n) == 0)
                     stop(kwdError, call.=FALSE)
                 description(object) <- value
                 return(object)
             })

setReplaceMethod("keyword", signature=c("flowFrame", "ANY"),
                 definition=function(object, value)
                 stop(kwdError, call.=FALSE))



## ==========================================================================
## plot method: We actually need to attach flowViz to do the plotting
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot",
          signature=signature(x="flowFrame",
                              y="ANY"),
          function(x, y, ...)
      {
          message("For plotting, please attach the 'flowViz' package.\n",
                  "   i.e., 'library(flowViz)'")
      })



## ==========================================================================
## nrow method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("nrow",
          signature=signature(x="flowFrame"),
          definition=function(x)
            return(nrow(x@exprs))
          )



## ==========================================================================
## ncol method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ncol",
          signature=signature(x="flowFrame"),
          definition=function(x)
            return(ncol(x@exprs))
          )



## ==========================================================================
## dim method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("dim",
          signature=signature(x="flowFrame"),
          definition=function(x)
      {
          d <- dim(x@exprs)
          names(d) <- c("events", "parameters")
          return(d)
      })



## ==========================================================================
## accessor method for feature names
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("featureNames",
          signature=signature(object="flowFrame"),
          definition=function(object)
          object@parameters$name
          )



## ==========================================================================
## accessor and replace methods for colnames of exprs slot
## this will also update the annotation in the parameters slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",
          signature=signature(x="flowFrame"),
          definition=function(x, do.NULL="missing", prefix="missing")
          colnames(exprs(x))
          )

setReplaceMethod("colnames",
                 signature=signature(x="flowFrame",
                                     value="ANY"),
                 definition=function(x, value)
             {
                 if(length(value) != ncol(x))
                     stop("colnames don't match dimensions of data matrix",
                          call.=FALSE)
                 pars <- parameters(x)
                 colnames(x@exprs) <- structure(as.character(value),
                                               names=rownames(pData(pars)))
                 pars$name <- value
                 names(pars$name) <- rownames(pData(pars))
                 x@parameters <- pars
                 return(x)
             })



## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
# setMethod("compensate",
#           signature=signature(x="flowFrame",
#                               spillover="matrix"),
#           definition=function(x, spillover, inv=TRUE, ...)
#       {
#           ## Make sure we're normalized to [0,1] and then invert
#           cols = colnames(spillover)
#           sel <- cols %in% colnames(x)
#           if(!all(sel))
#               stop("The following parameters in the spillover matrix\n are",
#                    " not present in the flowFrame:\n",
#                    paste(cols[!sel], collapse=", "), call.=FALSE)
#           e    = exprs(x)
#           if(inv)
#               e[,cols] = e[,cols]%*%solve(spillover/max(spillover))
#           else
#               e[,cols] = e[,cols]%*%spillover
#           exprs(x) = e
#           x
#       })
# 
# setMethod("compensate",
#           signature=signature(x="flowFrame",
#                               spillover="compensation"),
#           function(x, spillover, inv=TRUE, ...)
#       {
#           compensate(x, spillover=spillover@spillover,
#                      inv=spillover@invert)
#       })

setMethod("compensate",
          signature=signature(x="flowFrame",compensation="matrix"),
          definition=function(x,compensation,...)
      {    
          spillover=compensation
          cols = colnames(spillover)
          sel <- cols %in% colnames(x)
          if(!all(sel))
              stop("The following parameters in the spillover matrix\n are",
                   " not present in the flowFrame:\n",
                   paste(cols[!sel], collapse=", "), call.=FALSE)
          e = exprs(x)
          e[,cols]=t(solve(spillover)%*%t(e[,cols]))
          exprs(x) = e
          x
      })

setMethod("compensate",
          signature=signature(x="flowFrame",compensation="compensation"),
          function(x,compensation, ...)
          {   
              parameters=slot(compensation,"parameters")
              len=length(parameters)
              charParam=list()
              data=exprs(x)
    
              while(len>0)
              {
                  if(class(parameters[[len]])!="unitytransform")  ## process all transformed parameters
                  {
                      
                      if(class(parameters[[len]])=="transformReference")
                      { 
                          newCol=matrix(resolveTransformReference(parameters[[len]],data))
                          #newCol=eval(parameters[[len]])(data)
                          #newCol=matrix(eval(eval(parameters[[len]]))(data))
                          #newCol=matrix(eval(eval(parameters[[len]]))(data))
                          #colnames(newCol)=sprintf("_NEWCOL%03d_",len) 
                          colnames(newCol)=slot(parameters[[len]],"transformationId")
                          data <- cbind(data, newCol)
                      }
                      
                  } 
                  else 
                  {
                      charParam[[len]]=slot(parameters[[len]],"transformationId")
                  
                  }                
                  len=len-1                               
              }
            
              x<-flowFrame(data)
              spillover=compensation
              compensate(x,compensation@spillover,... )
          }
        )



## ==========================================================================
## transform method
## we are also making sure that the values of the dynamic range in the
## parameters slot are transformed accordingly. Note that this does not
## happen in the inline form of transformation!
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, ...)
      {
          e <- substitute(list(...))
          x <- `_data`
          par <- parameters(x)
          ranges <- range(x)
          tranges <- as.matrix(transform(as.data.frame(ranges),...))
          transformed <- as.matrix(transform(as.data.frame(exprs(x)),...))
          nc <- colnames(transformed)[-c(1:ncol(x@exprs))]
          colnames(transformed) <- c(colnames(x), nc)
          if(ncol(transformed) > ncol(x@exprs)) {
              ## Add any new parameter values, there might be a more elegant
              ## way to do that other than poking around in the deparsed
              ## arguments...
              args <- strsplit(gsub("^list\\(", "",
                                    unlist(strsplit(deparse(e), ","))), "=")
              nCol <- gsub(" *\"|\" *| *'|' *|` *|` *", "",
                           sapply(args,function(x) x[1]))
              nCol <- gsub(" ", "", nCol)
              oCol <- sapply(strsplit(sapply(args, function(x) x[2]),
                                      "`"), function(x) x[2])
              names(oCol) <- make.names(nCol)
              oparFrame <- pData(par)
              mt <- match(oCol[nc], oparFrame$name)
              nparFrame <- data.frame(name=names(oCol[nc]),
                                      desc=paste("derived from ",
                                      "transformation of", oCol[nc]),
                                      range=oparFrame[mt, "range"],
                                      minRange=tranges[1, nc],
                                      maxRange=tranges[2,nc])
              pData(par) <- rbind(oparFrame, nparFrame)
              }
              else {
                  cnames <- colnames(x)
                  par$minRange <- tranges[1,]
                  par$maxRange <- tranges[2,]
              }
              
              
              new("flowFrame",
                  exprs=transformed, 
                  parameters=par,#[,params$name],
                  description=description(x))
          })



## ==========================================================================
## filter method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## apply filter
setMethod("filter",
          signature=signature(x="flowFrame",
                              filter="filter"),
          definition=function(x, filter)
      {
          temp <- resolveTransforms(x, filter)
          x <- temp[["data"]]
          filter <- temp[["filter"]]
          allPar <- parameters(filter) %in% colnames(x)
          if(!all(allPar))
              stop("The following parameter(s) are not present in this ",
                   "flowFrame:\n", paste("\t", parameters(filter)[!allPar],
                                         collapse="\n"), call.=FALSE)
          result <- as(x %in% filter, "filterResult")
          identifier(result) <- identifier(filter)
          filterDetails(result, identifier(filter)) <- filter
          result@frameId <- identifier(x)
          result
      })

## apply filterSet
setMethod("filter",
          signature=signature(x="flowFrame",
                              filter="filterSet"),
          definition=function(x,filter)
      {
        ## A list of filter names sorted by dependency
          fl <- sort(filter,dependencies=TRUE)
          e <- filterSet()
          m <- as(filter,"list")
          r <-  lapply(fl,function(n) {
              e[[n]] <- eval(m[[n]])
              filter(x,e[n])
          })
          ## A place to evaluate filters
          ## e  = new.env()
          ## assign(".__frame",x,env=e)
          ## Evaluate each 
          ## r = lapply(as(filter,"list")[fl],function(f) {
          ## if(is(f,"formula")) f = eval(as.call(c(as.symbol("as"),
          ##      f,"filter")),e)
          ##     print(f)
          ## eval(as.call(c(as.symbol("filter"),as.symbol(".__frame"),f)),e)	
          ## })
          if(all(sapply(r,is,"logicalFilterResult"))) {
              ## Combine into a many filter result
              manyFilterResult(r,frameId=identifier(x),attr(fl,'AdjM'))
          } else 
          r
      })



## ===========================================================================
## spillover method
## Flow frames can have a SPILL keyword that is often used to store the
## spillover matrix. Attempt to extract it.
## ---------------------------------------------------------------------------
setMethod("spillover",
          signature=signature(x="flowFrame"),
          definition=function(x)
      {
          present <- keyword(x, c("spillover", "SPILL"))
          if(all(sapply(present, is.null)))
              stop("No spillover matrix stored in that flowFrame")
          else present
      })



## ==========================================================================
## wrappers for apply (over exprs slot)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("each_row",
          signature=signature(x="flowFrame"),
          definition=function(x,FUN,...)
      {
          apply(exprs(x),1,FUN,...)
      })

setMethod("each_col",
          signature=signature(x="flowFrame"),
          definition=function(x,FUN,...)
      {
          apply(exprs(x),2,FUN,...)
      })



## ==========================================================================
## Subset methods for flowFrame
## Why do we need the 'select' parameter? Wouldn't this be equivalent:
## Subset(x[,c(1,3)], subset)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
setMethod("Subset",
          signature=signature(x="flowFrame",
                              subset="filter"),
          definition=function(x,subset,select,...)
      {
          if(!missing(select))
              Subset(x,x %in% subset,select,...)
          else
              Subset(x,x %in% subset,...)
      })

setMethod("Subset",
          signature=signature(x="flowFrame",
                              subset="logical"),
          definition=function(x,subset,select,...)
      {
          if(!missing(select))
              x[subset & !is.na(subset),select]
          else x[subset,]
      })



## ==========================================================================
## range method for flowFrame
## get the range of possible data values from the parameters slot. Note that
## this range is not necessarily the range of the actual data values but
## it is the dynamic range of the measurement instrument
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Do we want the range to be set in a way to allow negative values?
## Need to fix that in read.FCS if yes
setMethod("range",
          signature=signature(x="flowFrame"),
          definition=function(x, ..., na.rm)
      {
          pind <- 1:nrow(parameters(x))
          channel <- unlist(list(...))
          if(length(channel)>0){
              if(!all(is(channel, "character")))
                  stop("All further arguments must be characters",
                       call.=FALSE)
                  pind <- match(channel, parameters(x)$name, nomatch=0)
          }
          if(any(pind==0))
              stop("'", paste(channel[which(pind==0)],
                              collapse="' and '", sep=""),
                   "' is/are no valid parameter(s) in this frame",
                   call.=FALSE) 
          ret <- rbind(min=parameters(x)$minRange[pind],
                       max=parameters(x)$maxRange[pind])
              colnames(ret) <- parameters(x)$name[pind]
          return(as.data.frame(ret))
      })



## ==========================================================================
## check for equality between two flowFrames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("==",
          signature=signature(e1="flowFrame",
                              e2="flowFrame"),
          definition=function(e1, e2)
      {
          return(all(as.numeric(exprs(e1)) == as.numeric(exprs(e2))) &&
                 all(as.character(pData(parameters(e1))) ==
                     as.character(pData(parameters(e2)))))
      })



## ==========================================================================
## show first or last values in the exprs matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("head",
          signature=signature(x="flowFrame"),
          definition=function(x, ...) head(exprs(x), ...))

setMethod("tail",
          signature=signature(x="flowFrame"),
          definition=function(x, ...) tail(exprs(x), ...))
          


## ==========================================================================
## comparison operators, these basically treat the flowFrame as a numeric
## matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("<",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) < e2)

setMethod(">",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) > e2)

setMethod("<=",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) <= e2)

setMethod(">=",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) >= e2)
