## ==========================================================================
## show and print methods display details about an object, and print usually
## allows for a bit mor fine-grained control.
## ==========================================================================






## ==========================================================================
## flowFrame
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
setMethod("show",
          signature=signature(object="filter"),
          definition=function(object)
          cat(paste("A filter named '", object@filterId, "'\n", sep="")))



## ==========================================================================
## filterReference
## ---------------------------------------------------------------------------
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
setMethod("show",
          signature=signature(object="transformFilter"),
          definition=function(object)
      {
          cat("transformed filter '", identifier(object), "'\n", sep="")
      })



## ==========================================================================
## transformMap
## ---------------------------------------------------------------------------
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
setMethod("show",
          signature=signature(object="filterResult"),
          definition=function(object)
          cat(paste("A filterResult produced by the filter named '",
                    object@filterId, "'\n", sep="")))



## ==========================================================================
## manyFilterResult
## --------------------------------------------------------------------------
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

setMethod("show",
		signature=signature(object="filtersList"),
		definition=function(object)
		{
			cat(sprintf("A list of %d filters .\n",length(object)))
		})

## ==========================================================================
## filterSet
## --------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="filterSet"),
          definition=function(object)
      {
          cat("A set of filter objects:\n")
          cat(paste(sapply(names(object), function(i)
                           identifier(object@env[[i]])), sep=",",
                    collapse=","))
          cat("\n")
      })



## ==========================================================================
## filterSummary
## --------------------------------------------------------------------------
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
setMethod("show",
          signature=signature(object="norm2Filter"),
          definition=function(object)
          {
            parms <- as.character(parameters(object))
            na  <-  is.na(parms)
            if(any(na))
              parms[na] <- "internal transformation"
            cat(ifelse(length(object@transformation), "transformed", ""),
                "norm2Filter '", identifier(object),
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
## workFlow. Recursively print all views in the workflow tree adding
## indentation
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="workFlow"),
          definition=function(object){
              cat("A flow cytometry workflow called '", object@name,
                  "'\n", sep="")
              tree <- get(object@tree)
              if(!is.null(tree)){
                nodes <- nodes(tree)
                if(length(nodes)>0){
                  traverseShow <- function(g, node, object, indent=1)
                  {
                      if(length(node)){
                          children <- unlist(adj(g, node))
                          for(i in children){
                              cat("\n")
                              print(get(i, object), indent=indent,
                                    parent=FALSE)
                              traverseShow(g, i, object, indent=indent+1)
                          }
                      }
                  }
                  cat("The following data views are provided:\n\n")
                  print(get(nodes[1], object))
                 traverseShow(tree, nodes[1], object)
                }
              }
              if(length(nodes)==0)
                  cat(" There is no view specified.\n")
          })



## ==========================================================================
## view. The print method allows more fine-grained control over the output,
## e.g., indentation
## ---------------------------------------------------------------------------
setMethod("print",
          signature=signature(x="view"),
          definition=function(x, indent=0, parent=TRUE)
      {
          b <- is.null(action(x))
          ind <- paste(rep("     ", indent), collapse="")
          if(indent>0)
              ind <- paste(ind, "  ", collapse="")
          cat(ind, ifelse(b, " Basic view", "View"), " '", x@name,
              "'\n", sep="")
          ## Only add the ID when the alias is non-unique
          if(!uniqueAlias(names(x), x))
              cat(ind, " (ID=", x@ID, ")\n", sep="")
          if(b){
              cat(ind, " on a ", class(Data(x)), " \n", sep="")
              cat(ind, " not associated to a particular ",
                  "action item\n", sep="")
          }else{
              cat(ind, " on a ", class(Data(parent(action(x)))),
                  " linked to \n", sep="")
              cat(ind, " ", sep="")
              print(action(x), indent=indent, parent=parent)
          }
      })

setMethod("show",
          signature=signature(object="view"),
          definition=function(object) print(object))




## ==========================================================================
## compensateActionItem. The print method allows more fine-grained control
## over the output, e.g., indentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("print",
          signature=signature(x="compensateActionItem"),
          definition=function(x, indent=0, parent=TRUE){
              ind <- paste(rep("     ", indent), collapse="")
              if(indent>0)
                  ind <- paste(ind, "  ", collapse="")
              cat("compensation action item '", x@name, "'\n", sep="")
              ## Only add the ID when the alias is non-unique
              if(!uniqueAlias(names(x), x))
                  cat(ind, " (ID=", x@ID, ")\n", sep="")
              if(parent)
                  cat(ind, " applied to view '", get(x@parentView)@name,
                      "' (ID=",identifier(x@parentView), ")\n", sep="")
              })

setMethod("show",
          signature=signature(object="compensateActionItem"),
          definition=function(object) print(object))



## ==========================================================================
## gateActionItem. The print method allows more fine-grained control
## over the output, e.g., indentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("print",
          signature=signature(x="gateActionItem"),
          definition=function(x, indent=0, parent=TRUE){
              ind <- paste(rep("     ", indent), collapse="")
              if(indent>0)
                  ind <- paste(ind, "  ", collapse="")
              cat("gate action item '", x@name, "'\n", sep="")
              ## Only add the ID when the alias is non-unique
              if(!uniqueAlias(names(x), x))
                  cat(ind, " (ID=", x@ID, ")\n", sep="")
              if(parent)
                  cat(ind, " applied to view '", get(x@parentView)@name,
                      "' (ID=",identifier(x@parentView), ")\n", sep="")
              })

setMethod("show",
          signature=signature(object="gateActionItem"),
          definition=function(object) print(object))



## ==========================================================================
## transformActionItem. The print method allows more fine-grained control
## over the output, e.g., indentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("print",
          signature=signature(x="transformActionItem"),
          definition=function(x, indent=0, parent=TRUE){
              ind <- paste(rep("     ", indent), collapse="")
              if(indent>0)
                  ind <- paste(ind, "  ", collapse="")
              cat("transform action item '", x@name, "'\n", sep="")
              ## Only add the ID when the alias is non-unique
              if(!uniqueAlias(names(x), x))
                  cat(ind, " (ID=", x@ID, ")\n", sep="")
              if(parent)
                  cat(ind, " applied to view '", get(x@parentView)@name,
                      "' (ID=",identifier(x@parentView), ")\n", sep="")
              })

setMethod("show",
          signature=signature(object="transformActionItem"),
          definition=function(object) print(object))





## ==========================================================================
## normalizeActionItem. The print method allows more fine-grained control
## over the output, e.g., indentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("print",
          signature=signature(x="normalizeActionItem"),
          definition=function(x, indent=0, parent=TRUE){
              ind <- paste(rep("     ", indent), collapse="")
              if(indent>0)
                  ind <- paste(ind, "  ", collapse="")
              cat("normalization action item '", x@name, "'\n", sep="")
              ## Only add the ID when the alias is non-unique
              if(!uniqueAlias(names(x), x))
                  cat(ind, " (ID=", x@ID, ")\n", sep="")
              if(parent)
                  cat(ind, " applied to view '", get(x@parentView)@name,
                      "' (ID=",identifier(x@parentView), ")\n", sep="")
              })

setMethod("show",
          signature=signature(object="normalizeActionItem"),
          definition=function(object) print(object))



## ==========================================================================
## subsettingActionItem. The print method allows more fine-grained control
## over the output, e.g., indentation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("print",
          signature=signature(x="subsettingActionItem"),
          definition=function(x, indent=0, parent=TRUE){
              ind <- paste(rep("     ", indent), collapse="")
              if(indent>0)
                  ind <- paste(ind, "  ", collapse="")
              cat("subsetting action item '", x@name, "'\n", sep="")
              ## Only add the ID when the alias is non-unique
              if(!uniqueAlias(names(x), x))
                  cat(ind, " (ID=", x@ID, ")\n", sep="")
              if(parent)
                  cat(ind, " applied to view '", get(x@parentView)@name,
                      "' (ID=",identifier(x@parentView), ")\n", sep="")
              })

setMethod("show",
          signature=signature(object="subsettingActionItem"),
          definition=function(object) print(object))


## ==========================================================================
## fcNullReference
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="fcNullReference"),
          definition=function(object) cat("NULL reference\n"))



## ==========================================================================
## fcReference
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="fcReference"),
          definition=function(object)
      {
          if(!resolved(object))
              cat(class(object) ," to unresolved object '",
                  object@ID, "'\n", sep="")
          else
              cat(class(object), " to ", class(get(object)),
                  " object '", object@ID, "'\n", sep="")
      })


## ==========================================================================
## transform
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="transform"),
          definition=function(object)
      {
          cat("transform object '", identifier(object), "'\n", sep="")
      })


setMethod("show",
          signature=signature(object="unitytransform"),
          definition=function(object)
      {
          cat("unitytransform on parameter '", parameters(object), "'\n",
              sep="")
      })



## ==========================================================================
## subsetting
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="subsetting"),
          definition=function(object)
      {
          cat("subsetting object '", identifier(object), "' to samples:\n",
              paste(object@indices, collapse=", "), "\n", sep="")
      })
