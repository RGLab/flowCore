## ==========================================================================
## flowFrames are the basic data structures for flow cytometry data
## ==========================================================================






## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## When subsetting columns (i.e., parameters) we also want to clean up the
## $PnX keywords
subsetKeywords <- function(x, j)
{
    kw <- description(x)
    par.id <- as.integer(gsub("^\\$P", "", rownames(pData(parameters(x)))))
    
    if(is.logical(j))
        j <- which(j)
    if(is.character(j))
        j <- match(j, colnames(x))
    
    to.del <- par.id[-j]
    pat <- paste(to.del, collapse = "|")
    pat <- paste0("^\\$P(", pat , ")[A-Z]$")
    sel <- grep(pat, names(kw))
    
    if(length(sel))
        x@description <- kw[-sel]
    return(x)
}


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
          switch(1+missing(i)+2*missing(j),
             {
                 x <- subsetKeywords(x, j)
                 exprs(x) <- exprs(x)[i, j, ..., drop=FALSE]
                 x
             },
             {
                 x <- subsetKeywords(x, j)
                 exprs(x) <- exprs(x)[ , j, ..., drop=FALSE]
                 x
             },
             {
                 exprs(x) <- exprs(x)[i,  , ..., drop=FALSE]
                 x
             },
             {
                 exprs(x) <- exprs(x)[ ,  , ..., drop=FALSE]
                 x
             } )
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
          definition=function(object){
              object@exprs
          })



setReplaceMethod("exprs",
                 signature=signature(object="flowFrame",
                                     value="matrix"),
                 definition=function(object, value)
             {
                 if(!is.numeric(value))
                     stop("replacement value for the 'exprs' slot ",
                          "of 'flowFrame' objects must be numeric",
                          "matrix", call.=FALSE)
                 mt <- match(colnames(value), colnames(object))
                 if(length(mt)==0)
                     stop("colnames missing in replacement value. ",
                          "Unable to match data columns", call.=FALSE)
                 if(any(is.na(mt)))
                     stop("the following parameters are not defined in the ",
                          "'flowFame' object:\n\t",
                          paste(colnames(value)[is.na(mt)], collapse=", "),
                          "\nunable to replace", call.=FALSE)
                 object@parameters <- object@parameters[mt,,drop=FALSE]
                 tmp <- object@exprs
                 
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
              
              kw <- keyword(object)
              if(hideInternal){
                  sel <- grep("^\\$", names(kw))
                  kw <- kw[-sel]
              }
              kw
            })

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

## select or combine keywords by function
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="function"),
          definition=function(object,keyword, ...)
          keyword(object, ...)
          )

## select keywords by combination of name and/or function
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="list"),
          definition=function(object,keyword)
      {
          sapply(keyword, function(k) {
              if(is.character(k))
                  object@description[k]
              else if(is.function(k))
                  k(object)
              else NA
          })
      })

## this is equivalent to the description method
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="missing"),
          function(object, compact = FALSE)
          {           
                       
                       desc <- object@description
                       if(compact)
                         desc <- kwfilter(desc)
                       desc
                       
          })
kwfilter <- function(desc){
  kn <- names(desc)
  pattern <- '(\\$)|(LASER)|(^P[1-9]{1,2})|(^FJ_)|(FCS)|FSC ASF|(CYTOMETER)|COMPENSATION|WINDOW|THRESHOLD|(CST )|SPILL|EXPORT |CREATOR|AUTOBS'
  kn <- kn[!grepl(pattern, kn)]  
  desc[kn]
}      

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
            return(nrow(exprs(x)))
          )



## ==========================================================================
## ncol method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ncol",
          signature=signature(x="flowFrame"),
          definition=function(x)
            return(ncol(exprs(x)))
          )



## ==========================================================================
## dim method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("dim",
          signature=signature(x="flowFrame"),
          definition=function(x)
      {
          d <- dim(exprs(x))
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
#' @rdname markernames
setGeneric("markernames",function(object,...) standardGeneric("markernames"))
#' get or update the marker names 
#' 
#' marker names corresponds to the 'desc' column of the phenoData of the flowFrame.
#' 
#' When extract marker names from a flowSet, it throws the warning if the marker names are not all the same across samples.
#' 
#' @param object flowFrame or flowSet
#' @param ... not used
#' @return marker names as a character vector. The marker names for FSC,SSC and Time channels are automatically excluded in the returned value.
#' When object is a flowSet and the marker names are not consistent across flowFrames, it returns a list of unique marker sets.
#' @rdname markernames
#' @examples 
#' 
#' data(GvHD)
#' fr <- GvHD[[1]]
#' markernames(fr)
#' 
#' chnls <- c("FL1-H", "FL3-H")
#' markers <- c("CD15", "CD14")
#' names(markers) <- chnls
#' markernames(fr) <- markers
#' markernames(fr)
#' 
#' fs <- GvHD[1:3]
#' markernames(fs)
#' @export
setMethod("markernames",
    signature=signature(object="flowFrame"),
    definition=function(object){
      param <- parameters(object)
      markers <- as.vector(param[["desc"]])
      ind <- grepl("time|fsc|ssc", param[["name"]], ignore.case = TRUE)
      markers <- markers[!ind]
      markers[!is.na(markers)]
    })
#' @rdname markernames
setGeneric("markernames<-",function(object, value) standardGeneric("markernames<-"))
#' @rdname markernames
#' @param value a named list or character vector. the names corresponds to the name(channel) and actual values are the desc(marker).
#' @export
setReplaceMethod("markernames",
                 signature=signature(object="flowFrame", value="ANY"), function(object, value){
                   
                   if(!is.character(value)){
                     stop("value must be a named character vector!")
                   }else{
                     chnls <- names(value)
                     if(is.null(chnls))
                       stop("value must be a named character vector!")
                     ind <- match(chnls, colnames(object))
                     misMatch <- is.na(ind)
                     if(any(misMatch))
                       stop("channel names not found in flow data: ", paste0(chnls[misMatch], collapse = ","))
                     
                     #convert factor to character  to avoid incorrect updating
                     old_desc <- pData(parameters(object))[,"desc"]
                     ids <- names(old_desc)[ind]
                     
                     if(is(old_desc, "factor"))
                     {
                       pData(parameters(object))[,"desc"] <- as.vector(old_desc)  
                     }
                     
                     pData(parameters(object))[ind, "desc"] <- as.vector(value)
                     
                     #restore to factor
                     if(is(old_desc, "factor")){
                       pData(parameters(object))[,"desc"] = factor(pData(parameters(object))[,"desc"])
                     }
                     
                     #update keyword
                     keyword(object)[ids] <- value
                     object
                   }
                   
                   
                 })

## ==========================================================================
## accessor and replace methods for colnames of exprs slot
## this will also update the annotation in the parameters slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("colnames",
          signature=signature(x="flowFrame"),
          definition=function(x, do.NULL="missing", prefix="missing")
          as.vector(parameters(x)[["name"]])
          )

setReplaceMethod("colnames",
                 signature=signature(x="flowFrame",
                                     value="ANY"),
                 definition=function(x, value)
             {
                 if(length(value) != ncol(x))
                     stop("colnames don't match dimensions of data matrix",
                          call.=FALSE)
                 value <- as.character(value)
                 
                 pars <- parameters(x)
                 ids <- rownames(pData(pars))
                 #update exprs
                 tmp <- exprs(x)
                 colnames(tmp) <- structure(value, names=ids)
                 #update param
                 pars$name <- value
                 names(pars$name) <- ids
                 x@parameters <- pars
                 exprs(x) <- tmp
                 #update keywords
                 
                 keys <- paste0(ids, "N")
                 keyword(x)[keys] <- value
                 return(x)
             })



## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
setMethod("compensate",
          signature=signature(x="flowFrame",
                              spillover="data.frame"),
          function(x, spillover)
      {
          compensate(x, as.matrix(spillover))
      })

setMethod("compensate",
          signature=signature(x="flowFrame",
                              spillover="matrix"),
          definition=function(x, spillover)
      {
#          checkClass(inv, "logical", 1)
#          if(inv)
#               spillover <- solve(spillover/max(spillover))
          cols <- colnames(spillover)
          sel <- cols %in% colnames(x)
          if(!all(sel))
              stop(keyword(x)[["FILENAME"]], "\nThe following parameters in the spillover matrix\n are",
                   " not present in the flowFrame:\n",
                   paste(cols[!sel], collapse=", "), call.=FALSE)
          e <- exprs(x)
		  
          e[, cols] <- e[, cols] %*% solve(spillover)
          exprs(x) = e
          x
      })

setMethod("compensate",
          signature=signature(x="flowFrame",
                              spillover="compensation"),
          function(x, spillover)
      {
          for(p in spillover@parameters)
              if(is(p, "transformReference"))
                  x  <- cbind2(x, matrix(resolveTransformReference(p, data),
                               dimnames=list(NULL, p@transformationId)))
          compensate(x, spillover@spillover)
      })




# transform method
# we are also making sure that the values of the dynamic range in the
# parameters slot are transformed accordingly. Note that this does not
# happen in the inline form of transformation!
# 
# @param _data flowFrame object
# @param translist transformList object. The recommended way of using this method.
# @param ... other arguments. e.g. `FL1-H` = myFunc(`FL1-H`)
#            but this form is not intended to be used in programmatic way since use non-standard evalution could fail to find
#            'myFunc' definition.   
setMethod("transform",
          signature=signature(`_data`="flowFrame"),
          definition=function(`_data`, translist, ...)
      {
        
          if(!(missing(translist))){
            #check if it is a transformList
            res <- try(class(translist), silent = TRUE)
            if(res != "transformList"){
              err_msg <- attr(res, "condition")[["message"]]
              err_msg <- paste(err_msg, "!Please make sure the unnamed argument is a valid 'transformList' object!")
              stop(err_msg)
            }else
                return(translist %on% `_data`)
          }else# dispach to .transform for named argument, assuming it is like `FSC-H`=asinhTrans(`FSC-H`) 
            .transform(`_data`, ...)
   
      })

#' take formal of transform(fs, `FSC-H`=asinhTrans(`FSC-H`))
#' which do the lazy evaluation
.transform <- function(`_data`, ...){
      e <- substitute(list(...))
      x <- `_data`
      par <- parameters(x)
      ranges <- range(x)
      tranges <- as.matrix(transform(as.data.frame(ranges),...))
      transformed <- as.matrix(transform(as.data.frame(exprs(x)),...))
      nc <- colnames(transformed)[-c(1:ncol(exprs(x)))]
      colnames(transformed) <- c(colnames(x), nc)
      if(ncol(transformed) > ncol(exprs(x))) {
        ## Add new parameter descriptions if there are any
        oCol <- lapply(e[2:length(e)], all.vars)
        nCol <- nc
        oparFrame <- pData(par)
        mt <- match(oCol[nc], oparFrame$name)
        nparFrame <- data.frame(name=nc,
            desc=paste("derived from",
                "transformation of",
                sapply(oCol[nc], paste,
                    collapse=" and ")),
            range=oparFrame[mt, "range"],
            minRange=tranges[1, nc],
            maxRange=tranges[2,nc])
        rownames(nparFrame) <- sprintf("$P%d", seq_len(nrow(nparFrame))+
                nrow(oparFrame))
        pData(par) <- rbind(oparFrame, nparFrame)
      }
      else
      {
        cnames <- colnames(x)
        par$minRange <- tranges[1,]
        par$maxRange <- tranges[2,]
      }
      #setting this flag is important 
      #since it tells the FCS parser that the latest range info needs to be read from flowCore_$P
      #instead of $PnR
      desc <- updateTransformKeywords(x)
      
      new("flowFrame", exprs = transformed, parameters = par,
          description = desc)
}

#' modify description to reflect the transformation
#' Involve inserting/updating 'transformation' and flowCore_$PnRmax keywords
#' @param fr flowFrame
#' @return updated description slot
updateTransformKeywords <- function(fr)
{
  description(fr) <- list(transformation="custom")
  desc <- description(fr)
  pd <- pData(parameters(fr))
  for(p in seq_along(pd[["name"]]))
  { 
    desc[[sprintf("flowCore_$P%sRmax", p)]] <- pd[p, "maxRange"]
    desc[[sprintf("flowCore_$P%sRmin", p)]] <- pd[p, "minRange"]
  }
  desc
}
## ==========================================================================
## filter method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## apply filter
setMethod("filter",
          signature=signature(x="flowFrame",
                              filter="filter"),
          definition=function(x, filter, method = "missing", sides = "missing", circular = "missing", init = "missing")
      {
            
          temp <- resolveTransforms(x, filter)
          x <- temp[["data"]]
          filter <- temp[["filter"]]
          allPar <- parameters(filter) %in% colnames(x)
          if (!is(filter, "setOperationFilter"))
		  {
              if(!all(allPar))
                  stop("The following parameter(s) are not present in this ",
                       "flowFrame:\n", paste("\t", parameters(filter)[!allPar],
                       collapse="\n"), call.=FALSE)
          }
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
          definition=function(x,filter, method = "missing", sides = "missing", circular = "missing", init = "missing")
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
          present <- keyword(x, .spillover_pattern)
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
        if (nrow(x) == 0) {
          x
        } else if(!missing(select)) {
          Subset(x, x %in% subset, select, ...)
        } else {
          Subset(x, x %in% subset, ...)
        }
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
          definition=function(x, ..., type = c("instrument", "data"), na.rm)
      {
          #deprecate the usage of dots
          args <- list(...)
          if(length(args)>0)
          {
            if(length(args)==1&&args[[1]]%in%c("instrument", "data"))
            {
              #overwrite type argument
              type <- args[[1]]
            }else
              stop("range method only accept two arguments: x(a flowFrame) and type = c('instrument', 'data')")
          }
            
            
          type <- match.arg(type)
          param <- parameters(x)
          pind <- 1:nrow(param)
          
          mir <- param$minRange[pind]
          mar <- param$maxRange[pind]
          
          if(is.null(mir)||is.null(mar)||type == "data")
            data_range <- apply(exprs(x), 2, range, na.rm=TRUE)
          
          if(is.null(mir)||type == "data")
              mir <- data_range[1,]
          
          if(is.null(mar)||type == "data")
              mar <- data_range[2,]
          
          ret <- rbind(min=mir, max=mar)
          colnames(ret) <- param$name[pind]
          
          
            
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



## ==========================================================================
## add data columns to a flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("cbind2",
          signature=signature(x="flowFrame", y="matrix"),
          definition=function(x, y)
      {
          nx <- nrow(x)
          ny <- nrow(y)
          if(mode(y) != "numeric")
              stop("'y' must be numeric matrix.", call.=FALSE)
          cn <- colnames(y)
          if(is.null(cn) || any(cn %in% colnames(x)))
              stop("Invalid or missing colnames in 'y'", call.=FALSE)
          exp <- cbind(exprs(x), y)
          parms <- parameters(x)
          parm <- pData(parms)
          range <- t(apply(y, 2, range))
          colnames(range) <- c("minRange", "maxRange")
          parm <- rbind(parm, data.frame(name=cn, desc=NA, range=NA, range))
          pData(parms) <- parm
          x@parameters <- parms
          exprs(x) <- exp
          for(i in seq_along(cn)){
              tmp <- list(cn[i])
              names(tmp) <- sprintf("$P%dN", i+ncol(x))
              description(x) <- tmp
          }
          return(x)
      })

setMethod("cbind2",
          signature=signature(x="flowFrame", y="numeric"),
          definition=function(x, y)
          stop("'y' has to be numeric matrix with colnames.", call.=FALSE))

## ==========================================================================
## get Index Sorted Data from a flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.getIndexSort<-function(x){
  if(class(x)!="flowFrame"){
    stop("x must be a flowFrame")
  }
  #extract keywords
  kw<-keyword(x)
  #subset keywords
  kw<-kw[grep("INDEX SORTING",names(kw))]
  if(length(kw)==0){
    message("No Index Sorting Information in this FCS")
    return(invisible(0))
  }
  .getCount<-function(x){
    cnt<-x[grep("COUNT",names(x))]
    lcnt<-length(cnt)
    try(cnt<-as.numeric(cnt[[1]]),silent=TRUE)
    if(is.na(cnt)||lcnt>1){
      stop("Invalid \"INDEX SORTING SORTED LOCATION COUNT\"")
    }
    cnt
  }
  .getDims<-function(x){
    dims<-x[grep("DIMENSION",names(x))]
    ldims<-length(dims)
    dims<-dims[[1]]
    if((ldims>1) || class(dims)!="character"){
      stop("Invalid \" INDEX SORTING DEVICE_DIMENSION\"")
    } 
    dims<-eval(parse(text=paste0("c(",dims,")")))
    return(dims)
  }
  .getLoc<-function(x){
    locs<-x[grep("LOCATIONS",names(x))]
    #locs<-locs[order(names(locs))]
	locids <- strsplit(names(locs),'_')
	locids <- unlist(lapply(locids,function(x) x[2]))
	locids <- as.numeric(locids)
	locs<-locs[order(locids)]
    m<-do.call(rbind,lapply(locs,function(y){
      m<-as.matrix(do.call(rbind,strsplit(strsplit(y,";")[[1]],",")))
      mode(m)<-"numeric"
      m
    }))
    return(m)
  }
  count<-.getCount(kw)
  dims<-.getDims(kw)
  locs<-.getLoc(kw)
  colnames(locs)<-c("XLoc","YLoc")
  mat<-cbind(exprs(x),locs)
  mat<-data.frame(name=keyword(x,"$FIL")[[1]],mat)
  return(mat)
}

setMethod("getIndexSort",signature="flowFrame",definition=function(x){
  .getIndexSort(x)
})

