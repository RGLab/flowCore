## ==========================================================================
## flowFrames are the basic data structures for flow cytometry data
## ==========================================================================



#' Append data columns to a flowFrame
#' 
#' Append data columns to a flowFrame
#' 
#' It is used to add extra data columns to the existing flowFrame.  It handles
#' keywords and parameters properly to ensure the new flowFrame can be written
#' as a valid FCS through the function \code{write.FCS} .
#' 
#' @name fr_append_cols
#' @aliases fr_append_cols
#' @usage 
#' fr_append_cols(fr, cols)
#' @param fr A \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param cols A numeric matrix containing the new data columns to be added.
#' Must has column names to be used as new channel names.
#' @return
#' 
#' A \code{\linkS4class{flowFrame}}
#' @author Mike Jiang
#' @keywords IO
#' @examples
#' 
#'   data(GvHD)
#'   tmp <- GvHD[[1]]
#'   
#'   kf <- kmeansFilter("FSC-H"=c("Pop1","Pop2","Pop3"), filterId="myKmFilter")
#'   fres <- filter(tmp, kf)
#'   cols <- as.integer(fres@subSet)
#'   cols <- matrix(cols, dimnames = list(NULL, "km"))
#'   tmp <- fr_append_cols(tmp, cols)
#'   
#'   tmpfile <- tempfile()
#'   write.FCS(tmp, tmpfile) 
#' 
#' 
#' @export
fr_append_cols <- function(fr, cols){
  checkClass(cols, "matrix")
  ncol <- ncol(cols)
  cn <- colnames(cols)
  if(length(cn) != ncol)
    stop("All columns in 'cols' must have colnames!")
  #add to pdata
  pd <- pData(parameters(fr))
  ncol_old <- ncol(fr)
  new_pid <- max(as.integer(gsub("\\$P", "", rownames(pd)))) + 1
  new_pid <- seq(new_pid, length.out = ncol)
  new_pid <- paste0("$P", new_pid)
  
  new_pd <- do.call(rbind, lapply(cn, function(i){
      vec <- cols[,i]
      rg <- range(vec)
      data.frame(name = i, desc = NA, range = diff(rg) + 1, minRange = rg[1], maxRange = rg[2])
    }))
  rownames(new_pd) <- new_pid
  pd <- rbind(pd, new_pd)  
  #add to exprs
  fr@exprs <- cbind(exprs(fr), cols)
  pData(parameters(fr)) <- pd
  #take care of flowCore_$PnRmax
  trans <- keyword(fr)[["transformation"]]
  if(!is.null(trans) && trans == "custom"){
    keyword(fr)[paste0("flowCore_", new_pid, "Rmax")] <- new_pd[new_pid, "maxRange"]
    keyword(fr)[paste0("flowCore_", new_pid, "Rmin")] <- new_pd[new_pid, "minRange"]
  }
  fr
}

## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## When subsetting columns (i.e., parameters) we also want to clean up the
## $PnX keywords
#' @examples 
#' fr <- GvHD[[1]]
#' fr <- subsetKeywords(fr, c(1,3,5))
#' kw <- keyword(fr)
#' kw[c("$P1N","$P2N", "$PAR")]
#' @noRd
subsetKeywords <- function(x, j)
{
    kw <- keyword(x)
    par.id <- as.integer(gsub("^\\$P", "", grep("^\\$P", rownames(pData(parameters(x))), value = TRUE)))
    if(is.logical(j))
      j <- which(j)
    if(is.character(j))
      j <- match(j, colnames(x))
    par.id <- par.id[j]
	  x@description <- filter_keywords(kw, par.id)
    return(x)
}

#' filter out $PnX keywords
#' 
#' @param kw a named list of keywords
#' @param par.id a vector of integer specifies the parameter ids to be perserved
#' @return a filtered list 
#' @examples 
#' data(GvHD)
#' fr <- GvHD[[1]]
#' kw <- keyword(fr)
#' kw <- filter_keywords(kw, c(1,3,5))
#' @export 
filter_keywords <- function(kw, par.id)
{
  par.id <- as.integer(par.id)
  stopifnot(is(par.id, "integer"))
  kn <- names(kw)
  
  all.pn <- kn[grep("^\\$P[0-9]+N$", kn)]
  all.pid <- as.integer(gsub("\\$P|N", "", all.pn))
  
  to.del.pid <- all.pid[!all.pid %in% par.id]
	pat <- paste(to.del.pid, collapse = "|")
	pat <- paste0("^\\$P(", pat , ")[A-Z]$")
	sel <- grep(pat, kn)
	
	if(length(sel))
		kw <- kw[-sel]
	return(kw)
}

## by indices or logical vectors
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
setMethod("exprs",
          signature=signature(object="flowFrame"),
          definition=function(object){
              object@exprs
          })



#' @export
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
#' @export
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
#' @export
setMethod("description",
          signature=signature(object="flowFrame"),
          function(object, hideInternal=FALSE){
			  .Deprecated("keyword")
              kw <- keyword(object)
              if(hideInternal){
                  sel <- grep("^\\$", names(kw))
                  kw <- kw[-sel]
              }
              kw
            })

## replace description entries
descError <- "Replacement value must be a named list."
#' @export
setReplaceMethod("description",
                 signature=signature(object="flowFrame",
                                     value="list"),
                 definition=function(object, value)
             {
				 .Deprecated("keyword<-")
                 n <- names(value)
                 if(length(n) == 0)
                     stop(descError, call.=FALSE)
                 object@description[n] <- value
                 return(object) })

#' @export
setReplaceMethod("description",
                 signature=signature(object="flowFrame",
                                     value="ANY"),
                 definition=function(object, value)
                 stop(descError, call.=FALSE))


## ==========================================================================
## accessor methods for individual items in the description slot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' Methods to retrieve keywords of a flowFrame
#' 
#' 
#' Accessor and replacement methods for items in the description slot (usually
#' read in from a FCS file header). It lists the \code{keywords} and its values
#' for a flowFrame specified by a character vector. Additional methods for
#' \code{function} and \code{lists} exists for more programmatic access to the
#' keywords.
#' 
#' The \code{keyword} methods allow access to the keywords stored in the FCS
#' files, either for a \code{flowFrame} or for a list of frames in a
#' \code{flowSet}. The most simple use case is to provide a character vector or
#' a list of character strings of keyword names. A more sophisticated version
#' is to provide a function which has to take one mandatory argument, the value
#' of this is the \code{flowFrame}. This can be used to query arbitrary
#' information from the \code{flowFrames} \code{description} slot or even the
#' raw data. The function has to return a single character string. The
#' \code{list} methods allow to combine functional and direct keyword access.
#' The replacement method takes a named character vector or a named list as
#' input. 
#' 
#' @name keyword-methods
#' @aliases keyword keyword-methods keyword,flowFrame,missing-method
#' keyword,flowFrame,function-method keyword,flowFrame,character-method
#' keyword,flowFrame,list-method keyword<- keyword<-,flowFrame,list-method
#' keyword<-,flowSet,list-method keyword<-,flowFrame,ANY-method
#' keyword<-,flowFrame,character-method keyword,flowSet,list-method
#' keyword,flowSet,ANY-method
#' @docType methods
#' 
#' @usage keyword(object, keyword, ...)
#' 
#' @param object Object of class \code{\link{flowFrame}}.
#' @param keyword Character vector or list of potential keywords or function.
#' If missing all keywords are returned.
#' @param ...  compact: logical scaler to indicate whether to hide all the
#' cytometer instrument and laser settings from keywords.
#' 
#' @section Methods: 
#' \describe{ 
#' \item{keyword(object = "flowFrame", keyword = "character")}{Return values for 
#' all keywords from the \code{description} slot in \code{object} that match 
#' the character vector \code{keyword}.}
#' 
#' \item{keyword(object = "flowFrame", keyword = "function")}{Apply the function in
#' \code{keyword} on the \code{\link{flowFrame}} \code{object}.  The function
#' needs to be able to cope with a single argument and it needs to return a
#' single character string. A typical use case is for instance to paste
#' together values from several different keywords or to compute some statistic
#' on the \code{flowFrame} and combine it with one or several other keywords. }
#' 
#' \item{keyword(object = "flowFrame", keyword = "list")}{Combine characters and
#' functions in a list to select keyword values.}
#' 
#' \item{keyword(object = "flowFrame", keyword = "missing")}{This is essentially an
#' alias for \code{\link{description}} and returns all keyword-value pairs.}
#' 
#' \item{keyword(object = "flowSet", keyword = "list")}{This is a wrapper around
#' \code{fsApply(object, keyword, keyword)} which essentially iterates over the
#' frames in the \code{\link{flowSet}}. }
#' 
#' \item{keyword(object = "flowSet", keyword = "ANY")}{This first coerces the
#' \code{keyword} (mostly a character vector) to a list and then calls the next
#' applicable method.  }
#' 
#' }
#' @author N LeMeur,F Hahne,B Ellis
#' @seealso \code{\link{description}}
#' @keywords methods
#' @examples
#' 
#' samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
#' keyword(samp)
#' keyword(samp, compact = TRUE)
#' 
#' keyword(samp, "FCSversion")
#' 
#' keyword(samp, function(x,...) paste(keyword(x, "SAMPLE ID"), keyword(x,
#' "GUID"), sep="_"))
#' 
#' keyword(samp)[["foo"]] <- "bar"
#' 
#' data(GvHD)
#' keyword(GvHD, list("GUID", cellnumber=function(x) nrow(x)))
#' 
#' keyword(GvHD)[["sample"]] <- sampleNames(GvHD))
#' 
#' 

## select keywords by name
#' @export
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="character"),
          function(object, keyword)
              structure(keyword(object, compact = FALSE)[keyword], names=keyword)
          )

## select or combine keywords by function
#' @export
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="function"),
          definition=function(object,keyword, ...)
          keyword(object, ...)
          )

## select keywords by combination of name and/or function
#' @export
setMethod("keyword",
          signature=signature(object="flowFrame",
                              keyword="list"),
          definition=function(object,keyword)
      {
          sapply(keyword, function(k) {
              if(is.character(k))
                  keyword(object, compact = FALSE)[k]
              else if(is.function(k))
                  k(object)
              else NA
          })
      })

## this is equivalent to the description method
#' @export
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
				 object@description <- value
				 return(object)
             })

setReplaceMethod("keyword", signature=c("flowFrame", "ANY"),
                 definition=function(object, value)
                 stop(kwdError, call.=FALSE))



## ==========================================================================
## plot method: We actually need to attach flowViz to do the plotting
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
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
#' @export
setMethod("nrow",
          signature=signature(x="flowFrame"),
          definition=function(x)
            return(nrow(exprs(x)))
          )



## ==========================================================================
## ncol method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("ncol",
          signature=signature(x="flowFrame"),
          definition=function(x)
            return(ncol(exprs(x)))
          )



## ==========================================================================
## dim method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
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
#' @export
setMethod("featureNames",
          signature=signature(object="flowFrame"),
          definition=function(object)
          object@parameters$name
          )
#' @rdname markernames
#' @export
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
      names(markers) <- param[["name"]]
      ind <- grepl("time|fsc|ssc", param[["name"]], ignore.case = TRUE)
      markers <- markers[!ind]
      markers[!is.na(markers)]
    })
#' @rdname markernames
#' @export
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
#' @noRd
#' @export
setMethod("colnames",
          signature=signature(x="flowFrame"),
          definition=function(x, do.NULL="missing", prefix="missing")
          as.vector(parameters(x)[["name"]])
          )
#' @export
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

## =====
## decompensate method
## =====

#' Decompensate a flowFrame
#' 
#' Reverse the application of a compensation matrix on a flowFrame
#'
#' @name decompensate
#' @docType methods
#' @aliases decompensate,flowFrame,matrix-method decompensate
#' decompensate,flowFrame,data.frame-method decompensate-methods
#' 
#' @usage decompensate(x, spillover)
#' @param x flowFrame. 
#' @param spillover matrix or data.frame. 
#'
#' @return a decompensated flowFrame
#'
#' @examples
#' library(flowCore)
#' f = list.files(system.file("extdata",
#'    "compdata",
#'    "data",
#'    package="flowCore"),
#'  full.name=TRUE)[1]
#' f = read.FCS(f)
#' spill = read.csv(system.file("extdata",
#'        "compdata", "compmatrix",
#'         package="flowCore"), 
#'         ,sep="\t",skip=2)
#' colnames(spill) = gsub("\\.","-",colnames(spill))
#' f.comp = compensate(f,spill)
#' f.decomp = decompensate(f.comp,as.matrix(spill))
#' sum(abs(f@exprs-f.decomp@exprs))
#' all.equal(decompensate(f.comp,spill)@exprs,decompensate(f.comp,as.matrix(spill))@exprs)
#' all.equal(f@exprs,decompensate(f.comp,spill)@exprs)
#'
#' @export
setMethod("decompensate",
			 signature = signature(x="flowFrame", spillover="matrix"),
			 function(x, spillover) 
			 {
		    cols <- .check_spillover(spillover, x)
			 	e <- exprs(x)
			 	e[, cols] <- e[, cols] %*% spillover
			 	exprs(x) = e
			 	x
			 }
)

#' @param spillover 
#' @param x flowFrame
#' @noRd
.check_spillover <- function(spillover, x){
  cols <- colnames(spillover)
  if(is.null(cols))
    stop("The spillover matrix must have colnames!")
  sel <- cols %in% colnames(x)
  if(!all(sel))
    stop(keyword(x)[["FILENAME"]], "\nThe following parameters in the spillover matrix\n are",
         " not present in the flowFrame:\n",
         paste(cols[!sel], collapse=", "), call.=FALSE)
  cols
  
}

#' @export
setMethod("decompensate",
			 signature = signature(x="flowFrame", spillover="data.frame"),
			 function(x, spillover) 
			 {
			 	spillover = as.matrix(spillover)
			 	decompensate(x = x, spillover = spillover)
			 }
)

## ===========================================================================
## compensate method
## ---------------------------------------------------------------------------
#' Compensate a flowFrame
#'
#' @name compensate
#' @aliases  compensate
#' @param x flowFrame. 
#' @param spillover matrix or data.frame. 
#'
#' @return a compensated flowFrame
#'
#' @export
#' @noRd
setMethod("compensate",
          signature=signature(x="flowFrame",
                              spillover="data.frame"),
          function(x, spillover)
      {
          compensate(x, as.matrix(spillover))
      })

#' @export
setMethod("compensate",
          signature=signature(x="flowFrame",
                              spillover="matrix"),
          definition=function(x, spillover)
      {
#          checkClass(inv, "logical", 1)
#          if(inv)
#               spillover <- solve(spillover/max(spillover))
          cols <- .check_spillover(spillover, x)            
          e <- exprs(x)
		  
          e[, cols] <- e[, cols] %*% solve(spillover)
          exprs(x) = e
          x
      })

#' @export
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




#' Transform a flowFrame or flowSet
#' 
#' Similar to the base transform method, this will transform the values of
#' a flowFrame or flowSet object according to the transformations specified
#' in one of two ways:
#' 1. a [transformList][flowCore::transformList-class] or list of [transform][flowCore::transform-class] objects
#' 2. named arguments specifying transformations to be applied to channels (see details)
#' 
#' @name transform
#' @md
#' @usage transform(_data, trans, ...)
#' @param _data a flowFrame or flowSet object
#' @param trans a transformList object
#' @param ... other arguments. e.g. `FL1-H` = myFunc(`FL1-H`)
#' 
#' @details To specify the transformations in the second way, the names of these arguments 
#' should correspond to the new channel names and the values should be functions applied to
#' channels currently present in the flowFrame or flowSet. There are a few examples below.
#' 
#' @examples 
#' data(GvHD)
#' # logarithmically transform FL1-H and FL2-H for the entire flowSet
#' # using a transformList
#' fs <- transform(GvHD,
#'                  transformList(c("FL1-H", "FL2-H"), list(log, log)))
#'                  
#' # transform a single flowFrame using named arguments. Note the first
#' # transformation will overwrite FL1-H while the second will create a new
#' # channel 
#' fr <- transform(GvHD[[1]],
#'                  `FL1-H`=log(`FL1-H`),
#'                  `logFL2`=log(`FL2-H`))
#' 
#' 
#' @export
# 
# We are also making sure that the values of the dynamic range in the
# parameters slot are transformed accordingly. Note that this does not
# happen in the inline form of transformation!
# 
# The `FL1-H` = myFunc(`FL1-H` ) form is not intended to be used in programmatic way
# since use non-standard evalution could fail to find 'myFunc' definition. 
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
#' @noRd
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
  keyword(fr)[["transformation"]] <- "custom"
  desc <- keyword(fr)
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

#' Filter FCS files
#' 
#' 
#' These methods link filter descriptions to a particular set of flow cytometry
#' data allowing for the lightweight calculation of summary statistics common
#' to flow cytometry analysis.
#' 
#' 
#' The \code{filter} method conceptually links a filter description,
#' represented by a \code{\linkS4class{filter}} object, to a particular
#' \code{\linkS4class{flowFrame}}. This is accomplished via the
#' \code{\linkS4class{filterResult}} object, which tracks the linked frame as
#' well as caching the results of the filtering operation itself, allowing for
#' fast calculation of certain summary statistics such as the percentage of
#' events accepted by the \code{filter}. This method exists chiefly to allow
#' the calculation of these statistics without the need to first
#' \code{\link{Subset}} a \code{\linkS4class{flowFrame}}, which can be quite
#' large.
#' 
#' When applying on a \code{flowSet}, the \code{filter} argument can either be
#' a single \code{filter} object, in which case it is recycled for all frames
#' in the set, or a named list of \code{filter} objects. The names are supposed
#' to match the frame identifiers (i.e., the output of \code{sampleNames(x)} of
#' the \code{flowSet}. If some frames identifiers are missing, the particular
#' frames are skipped during filtering. Accordingly, all \code{filters} in the
#' filter list that can't be mapped to the \code{flowSet} are ignored. Note
#' that all \code{filter} objects in the list must be of the same type, e.g.
#' \code{rectangleGates}.
#' 
#' @name filter-methods
#' @docType methods
#' @aliases filter filter-method filter,flowFrame-method filter,flowFrame,filter-method
#' filter,flowFrame,rectangleGate filter,flowFrame,polygonGate
#' filter,flowFrame,norm2Filter filter,flowSet,filter-method
#' filter,flowSet,list-method filter,flowSet,filterList-method
#' summary,filter-method show,filter-method length,filter-method
#' formula,filter-method character,filter-method name,filter-method
#' call,filter-method identifier<-,filter,character-method
#' 
#' @usage
#' filter(x, filter, method = c("convolution", "recursive"), 
#' sides = 2L, circular = FALSE, init = NULL)
#' 
#' @param x Object of class \code{\linkS4class{flowFrame}} or
#' \code{\linkS4class{flowSet}}.
#' @param filter An object of class \code{\linkS4class{filter}} or a named list
#' \code{filters}.
#' @param method,sides,circular,init These arguments are not used.
#' @return
#' 
#' A \code{\linkS4class{filterResult}} object or a
#' \code{\linkS4class{filterResultList}} object if \code{x} is a
#' \code{\linkS4class{flowSet}}. Note that \code{\linkS4class{filterResult}}
#' objects are themselves filters, allowing them to be used in filter
#' expressions or \code{Subset} operations.
#' @author F Hahne, B. Ellis, N. Le Meur
#' @seealso \code{\link{Subset}}, \code{\linkS4class{filter}}, \code{\linkS4class{filterResult}}
#' @keywords methods
#' @examples
#' 
#' ## Filtering a flowFrame
#' samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
#' rectGate <- rectangleGate(filterId="nonDebris","FSC-H"=c(200,Inf))
#' fr <- filter(samp,rectGate)
#' class(fr)
#' summary(fr)
#' 
#' ## filtering a flowSet
#' data(GvHD)
#' foo <- GvHD[1:3]
#' fr2 <- filter(foo, rectGate)
#' class(fr2)
#' summary(fr2)
#' 
#' ## filtering a flowSet using different filters for each frame
#' rg2 <- rectangleGate(filterId="nonDebris","FSC-H"=c(300,Inf))
#' rg3 <- rectangleGate(filterId="nonDebris","FSC-H"=c(400,Inf))
#' flist <- list(rectGate, rg2, rg3)
#' names(flist) <- sampleNames(foo)
#' fr3 <- filter(foo, flist)
#' 

## apply filter
#' @export
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


## ===========================================================================
## spillover method
## Flow frames can have a SPILL keyword that is often used to store the
## spillover matrix. Attempt to extract it.
## ---------------------------------------------------------------------------
#' @export
setMethod("spillover",
          signature=signature(x="flowFrame"),
          definition=function(x)
      {
          present <- keyword(x, .spillover_pattern)
          if(all(sapply(present, is.null)))
              stop("No spillover matrix stored in that sample")
          else present
      })



## ==========================================================================
## wrappers for apply (over exprs slot)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("each_row",
          signature=signature(x="flowFrame"),
          definition=function(x,FUN,...)
      {
          apply(exprs(x),1,FUN,...)
      })

#' @export
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
#' Subset a flowFrame or a flowSet
#' 
#' An equivalent of a \code{\link{subset}} function for
#' \code{\linkS4class{flowFrame}} or a \code{\linkS4class{flowSet}} object.
#' Alternatively, the regular subsetting operators can be used for most of the
#' topics documented here.
#' 
#' 
#' The \code{Subset} method is the recommended method for obtaining a
#' \code{\linkS4class{flowFrame}} that only contains events consistent with a
#' particular filter. It is functionally equivalent to
#' \code{frame[as(filter(frame,subset),"logical"),]} when used in the
#' \code{\linkS4class{flowFrame}} context. Used in the
#' \code{\linkS4class{flowSet}} context, it is equivalent to using
#' \code{\link{fsApply}} to apply the filtering operation to each
#' \code{\linkS4class{flowFrame}}.
#' 
#' Additionally, using \code{Subset} on a \code{\linkS4class{flowSet}} can also
#' take a named \code{list} as the subset. In this case, the names of the list
#' object should correspond to the \code{sampleNames} of the flowSet, allowing
#' a different filter to be applied to each frame. If not all of the names are
#' used or excess names are present, a warning will be generated but the valid
#' filters will be applied for the rare instances where this is the intended
#' operation. Note that a \code{\link{filter}} operation will generate a list
#' of \code{\linkS4class{filterResult}} objects that can be used directly with
#' \code{Subset} in this manner.
#' 
#' @name Subset-methods
#' @aliases Subset Subset,flowFrame,filter-method Subset,flowSet,ANY-method
#' Subset,flowFrame-method Subset,flowFrame,logical-method
#' Subset,flowSet,list-method Subset,flowSet,filterResultList-method
#' Subset,flowSet,ANY
#' 
#' @usage Subset(x, subset, ...)
#' 
#' @param x The flow object, frame or set, to subset.
#' @param subset A filter object or, in the case of \code{flowSet} subsetting,
#' a named list of filters.
#' @param \dots Like the original \code{\link{subset}} function, you can also
#' select columns.
#' 
#' @return
#' 
#' Depending on the original context, either a \code{\linkS4class{flowFrame}}
#' or a \code{\linkS4class{flowSet}}.
#' @author B. Ellis
#' @seealso \code{\link{split}}, \code{\link{subset}}
#' @keywords manip
#' @examples
#' 
#' sample <- read.flowSet(path=system.file("extdata", package="flowCore"),
#' pattern="0877408774")
#' result <- filter(sample, rectangleGate("FSC-H"=c(-Inf, 1024)))
#' result
#' Subset(sample,result)
#' 
#' 
#' @export
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

#' @export
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
#' @export
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
          if(type == "data")
            pnames <- colnames(exprs(x))[pind]
          else
            pnames <- param$name[pind]
          colnames(ret) <- pnames
          
          
            
          return(as.data.frame(ret))
      })



## ==========================================================================
## check for equality between two flowFrames
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
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
#' @export
setMethod("head",
          signature=signature(x="flowFrame"),
          definition=function(x, ...) head(exprs(x), ...))

#' @export
setMethod("tail",
          signature=signature(x="flowFrame"),
          definition=function(x, ...) tail(exprs(x), ...))
          


## ==========================================================================
## comparison operators, these basically treat the flowFrame as a numeric
## matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("<",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) < e2)

#' @export
setMethod(">",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) > e2)

#' @export
setMethod("<=",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) <= e2)

#' @export
setMethod(">=",
          signature=signature(e1="flowFrame",
                              e2="ANY"),
          definition=function(e1, e2) exprs(e1) >= e2)



## ==========================================================================
## add data columns to a flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("cbind2",
          signature=signature(x="flowFrame", y="matrix"),
          definition=function(x, y)
      {
		  .Deprecated("fr_append_cols")
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
			  kn <- sprintf("$P%dN", i+ncol(x))
              keyword(x)[kn] <- tmp
          }
          return(x)
      })

#' @export
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

#' Extract Index Sorted Data from an FCS File
#' 
#' 
#' Retrieve a data frame of index sorted data and sort indices from an FCS
#' file.
#' 
#' The input FCS file should already be compensated.  Index sorting permits
#' association of cell-level fluorescence intensities with downstream data
#' collection on the sorted cells. Cells are sorted into a plate with
#' \code{X,Y} coordinates, and those coordinates are stored in the FCS file.
#' 
#' This function will extract the data frame of flow data and the \code{X,Y}
#' coordinates for the cell-level data, which can be written to a text file, or
#' concatenated with sample-level information and analyzed in R. The
#' coordinates are names 'XLoc','YLoc', and a 'name' column is also prepended
#' with the FCS file name.
#' 
#' @name getIndexSort
#' @aliases getIndexSort getIndexSort-methods getIndexSort,flowFrame-method
#' @docType methods
#' @usage NULL
#' @return
#' 
#' Matrix of fluorescence intensities and sort indices for plate location.
#' When no index sorting data is available, invisibly returns 0. Test for 0 to
#' check success.
#' @section Methods:
#' 
#' \describe{
#' 
#' \item{getIndexSort(x = "flowFrame")}{Return a matrix of fluorescence intensities and
#' indices into the sorting plate for each cell.}
#' 
#' }
#' @author G. Finak
#' @keywords methods
#' @examples
#' 
#'  samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
#'  # This will return a message that no index sorting data is available
#'  getIndexSort(samp)
#' 
#' 
#' @export
setMethod("getIndexSort",signature="flowFrame",definition=function(x){
  .getIndexSort(x)
})

