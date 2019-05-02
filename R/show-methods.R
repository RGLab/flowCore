## ==========================================================================
## show and print methods display details about an object, and print usually
## allows for a bit mor fine-grained control.
## ==========================================================================






## ==========================================================================
## flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="flowFrame"),
          definition=function(object)
      {
          
          cat(paste(class(object), " object '", identifier(object),
                    "'\nwith ", nrow(object), " cells and ",
                    ncol(object), " observables:\n", sep=""))
          show(pData(parameters(object)))
          cat(paste(length(description(object)), " keywords are stored in the ",
                    "'description' slot\n", sep = ""))
          return(invisible(NULL))
      })



## ==========================================================================
## flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="flowSet"),
          definition=function(object)
      {
          cat("A flowSet with",length(object),"experiments.\n\n")
          
          pd <- phenoData(object)
          varDesc <- varMetadata(pd)$labelDescription
          varDesc <- varDesc[!is.na(varDesc)]
          
          if(length(varDesc) > 0){
            if(any(varDesc != "Name")){
              show(phenoData(object))
              cat("\n")
          }
        }
          cat("  column names:\n  ")
          cat(paste(object@colnames,sep=","))
          cat("\n")
      })



## ==========================================================================
## Compensation object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="compensation"),
          definition=function(object){
              cat("Compensation object '", object@compensationId,
                  "':\n", sep="")
              if(ncol(object@spillover)){
				  print(object@spillover)
#                   if(!object@invert)
#                       print(object@spillover)
#                   else
#                       print(solve(object@spillover/max(object@spillover)))
#                      ;
              }else{
                  cat("The spillover matrix is empty\n")
              }
          })



## ==========================================================================
## filter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="filter"),
          definition=function(object)
          cat(paste("A filter named '", object@filterId, "'\n", sep="")))



## ==========================================================================
## filterReference
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="filterReference"),
          definition=function(object)
      {
          if(exists(object@name,envir=object@env))
              cat(paste("A reference to a filter named '",
                        identifier(object), "'\n", sep=""))
          else
              cat(paste("An unresolvable reference to a filter named '",
                        object@name, "'\n", sep=""))
      })



## ==========================================================================
## complementFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="complementFilter"),
          definition=function(object)
      {
          cat("filter '", identifier(object),
              "', the complement of\n", sep="")
          print(object@filters[[1]])

      })



## ==========================================================================
## subsetFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="subsetFilter"),
          definition=function(object)
      {
          cat("filter '", identifier(object),
              "'\nthe filtering operation defined by\n", sep="")
          print(object@filters[[1]])
          cat("after subsetting by\n")
          print(object@filters[[2]])
      })


## ==========================================================================
## unionFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="unionFilter"),
          definition=function(object)
      {
          cat("filter '", identifier(object),
              "'\nthe union of the ", length(object@filters),
              " filters\n\n", sep="")
          for(i in 1:length(object@filters)){
              print(object@filters[[i]])
              cat("\n")
          }
      })



## ==========================================================================
## intersectFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="intersectFilter"),
          definition=function(object)
      {
          cat("filter '", identifier(object),
              "'\nthe intersection between the ", length(object@filters),
              " filters\n\n", sep="")
          for(i in 1:length(object@filters)){
              print(object@filters[[i]])
              cat("\n")
          }
      })



## ==========================================================================
## transformFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="transformFilter"),
          definition=function(object)
      {
          cat("transformed filter '", identifier(object), "'\n", sep="")
      })



## ==========================================================================
## transformMap
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature(object="transformMap"),
          definition=function(object)
      {
          cat("transformMap for parameter '",
              object@input, "' mapping to '",
              object@output, "'\n", sep="")
      })



## ==========================================================================
## filterResult
## --------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="filterResult"),
          definition=function(object)
          cat(paste("A filterResult produced by the filter named '",
                    object@filterId, "'\n", sep="")))



## ==========================================================================
## manyFilterResult
## --------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="manyFilterResult"),
          definition=function(object)
      {
          n <- names(object)
          cat("A filter result containing potentially overlapping populations:\n")
          cat(paste(n,collapse=", "))
          cat("\n")
      })



## ==========================================================================
## multipleFilterResult
## --------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="multipleFilterResult"),
          definition=function(object)
      {
          cat(paste("A filterResult produced by the filter named '",
                    object@filterId, "'\n resulting in multiple ",
                    "populations:\n", paste("\t", names(object), collapse="\n"),
                    "\n", sep=""))
      })



## ==========================================================================
## filterResultList
## --------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="filterResultList"),
          definition=function(object)
      {
          cat(paste("A list of filterResults for a flowSet containing",
                    length(object), "frames\nproduced by",
                    ifelse(length(object@filterId)>1,
                           "frame-specific filters\n",
                           paste("the filter named '",
                                 object@filterId, "'\n",
                                 sep=""))))
      })


## ==========================================================================
## filterList
## --------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="filterList"),
          definition=function(object)
      {
          cat(sprintf("A list of %d filters with filterId '%s'.\n",
                      length(object), identifier(object)))
      })

## ==========================================================================
## filters
## --------------------------------------------------------------------------
#' @export

setMethod("show",
		signature=signature(object="filters"),
		definition=function(object)
		{
			cat(sprintf("A list of %d filters applied to a flowFrame.\n",
							length(object)))
		})
## ==========================================================================
## filtersList
## --------------------------------------------------------------------------

#' @export

setMethod("show",
		signature=signature(object="filtersList"),
		definition=function(object)
		{
			cat(sprintf("A list of %d filters .\n",length(object)))
		})



## ==========================================================================
## filterSummary
## --------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="filterSummary"),
          definition=function(object)
      {
          for(i in seq(along=object@name)) {
              cat(sprintf("%s: %d of %d events (%.2f%%)\n",
                          object@name[i],
                          object@true[i], object@count,
                          object@p[i]*100))
          }
      })

## A bit mor control over the output (identation)
#' @export

setMethod("print",
          signature=signature(x="filterSummary"),
          definition=function(x, indent=0, verbose=TRUE)
      {
          if(verbose){
              for(i in seq(along=x@name)) {
                  cat(rep(" ", indent),
                      sprintf("%s: %d of %d events (%.2f%%)\n",
                              x@name[i],
                              x@true[i], x@count,
                              x@p[i]*100), sep="")
              }
          }
          return(invisible(x))
      })

## ==========================================================================
## expressionFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="expressionFilter"),
          definition=function(object)
      {
          msg <- paste("expression filter '", identifier(object),
                       "' evaluating the expression:\n",
                       paste(object@deparse, collapse="\n"), sep="")
          cat(msg)
          cat("\n")
          invisible(msg)
      })



## ==========================================================================
## ellipsoidGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="ellipsoidGate"),
          definition=function(object)
      {
        parms <- as.character(parameters(object))
        na  <-  is.na(parms)
        if(any(na))
          parms[na] <- "internal transformation"
        cat("Ellipsoid gate '", identifier(object),
            "' in dimensions ", sep="")
        cat(paste(parms, sep="", collapse=" and "))
        cat("\n")
      })



## ==========================================================================
## kmeansFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="kmeansFilter"),
          definition=function(object)
      {
        parms <- as.character(parameters(object))
        na  <-  is.na(parms)
        if(any(na))
          parms[na] <- "internal transformation"
        msg <- paste("k-means filter '", object@filterId,
                     "' in dimension ", parms[1],
                     "\nwith ", length(object), " populations (",
                     paste(object@populations, collapse=","),
                     ")", sep="")
        cat(msg)
        cat("\n")
        invisible(msg)
      })



## ==========================================================================
## norm2Filter
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="norm2Filter"),
          definition=function(object)
          {
            parms <- as.character(parameters(object))
            na  <-  is.na(parms)
            if(any(na))
              parms[na] <- "internal transformation"
            cat("norm2Filter '", identifier(object),
                "' in dimensions ", sep="")
            cat(paste(parms, sep="", collapse=" and "),
                "with parameters:\n")
            cat("  method:", object@method, "\n")
            cat("  scale.factor:", object@scale.factor, "\n")
            cat("  n:", object@n, "\n")
            cat("\n")
          })



## ==========================================================================
## polygonGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="polygonGate"),
          definition=function(object)
          {
            parms <- as.character(parameters(object))
            na  <-  is.na(parms)
            if(any(na))
              parms[na] <- "internal transformation"
            nb <-  nrow(object@boundaries)
            cat("Polygonal gate '", identifier(object) ,"' with ",
                ifelse(all(is.na(object@boundaries)), 0, nb),
                " vertices in dimensions ", sep="")
            cat(paste(parms, sep="", collapse=" and "))
            cat("\n")
          })



## ==========================================================================
## quadGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="quadGate"),
          definition=function(object)
      {
        parms <- as.character(parameters(object))
        na  <-  is.na(parms)
        if(any(na))
          parms[na] <- "internal transformation"
        cat("Quadrant gate '", identifier(object),
            "' with dimensions:\n", sep="")
        for(i in seq(along=parameters(object))) {
          cat("  ")
          cat(parms[i])
          cat(": ")
          cat(object@boundary[i])
          cat("\n")
        }
      })



## ==========================================================================
## rectangleGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export

setMethod("show",
          signature=signature(object="rectangleGate"),
          definition=function(object)
      {
          parms <- as.character(parameters(object))
          na  <-  is.na(parms)
          if(any(na))
            parms[na] <- "internal transformation"
          cat("Rectangular gate '", identifier(object),
              "' with dimensions:\n", sep="")
          for(i in seq_along(parms)){
              cat("  ")
              if(!is.na(parms[i]))
                  cat(parms[i])
              else
                  cat("anonymous parameter")
              cat(": (")
              cat(paste(object@min[i],object@max[i],sep=","))
              cat(")\n")
	}
      })



## ==========================================================================
## sampleFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="sampleFilter"),
          function(object)
      {
          msg <- paste("sample filter '", object@filterId,
                       "' returning objects with ", object@size," rows",
                       sep="")
          cat(msg)
          cat("\n")
          invisible(msg)
      })



## ==========================================================================
## timeFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="timeFilter"),
          definition=function(object)
      {
          msg <- paste("time filter '",object@filterId,
                       "' with settings:\n  bandwidth=",
                       object@bandwidth, sep="")
          cat(msg)
          if(length(object@binSize))
              cat("\n  binSize=", object@binSize, sep="")
          if(length(object@timeParameter))
              cat("\n  timeParameter=", object@timeParameter, sep="")
          cat("\n")
          invisible(msg)
      })


## ==========================================================================
## boundaryFilter
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="boundaryFilter"),
          definition=function(object)
      {
          msg <- paste("boundaryFilter '",object@filterId,
                       "' operating on ",
                       sprintf("channel%s\n", ifelse(length(object@side)==1, "", "s:")),
                       paste(" ", parameters(object), " (tolerance=",
                             signif(object@tolerance,3),
                             ", boundary=", object@side, ")\n", sep="",
                             collapse=""), sep="")
          cat(msg)
          invisible(msg)
      })



## ==========================================================================
## transform
## ---------------------------------------------------------------------------
#' @export

setMethod("show",
          signature=signature(object="transform"),
          definition=function(object)
      {
          cat("transform object '", identifier(object), "'\n", sep="")
      })


#' @export

setMethod("show",
          signature=signature(object="unitytransform"),
          definition=function(object)
      {
          cat("unitytransform on parameter '", parameters(object), "'\n",
              sep="")
      })
