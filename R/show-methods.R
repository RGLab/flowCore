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
          dm <- dim(exprs(object))
          cat(paste("flowFrame object '", identifier(object),
                    "'\nwith ", dm[1], " cells and ", 
                    dm[2], " observables:\n", sep=""))
          show(pData(parameters(object)))
          cat(paste(length(description(object)), " keywords are stored in the ",
                    "'descripton' slot\n", sep = ""))
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
          if(any(varMetadata(phenoData(object))$labelDescription != "Name")){
              show(phenoData(object))
              cat("\n")
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
                  if(!object@invert)
                      print(object@spillover)
                  else
                      print(solve(object@spillover/max(object@spillover)))
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
          if(exists(object@name,env=object@env))
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
          cat("A filter result containing: ")
          cat(paste(n,sep=","))
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
          definition=function(x, indent=0)
      {
          for(i in seq(along=x@name)) {
              cat(rep(" ", indent),
                  sprintf("%s: %d of %d events (%.2f%%)\n",
                          x@name[i],
                          x@true[i], x@count,
                          x@p[i]*100), sep="")
          }
          return(invisible(x))
      })



## ==========================================================================
## curv1Filter
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="curv1Filter"),
          definition=function(object)
      {
          msg <- paste("1D curvature filter '",object@filterId,
                       "' in dimension ",
                       object@parameters, "\nwith settings:",
                       "\n  bwFac=", object@bwFac, "\n  gridsize=",
                       paste(object@gridsize, collapse=",", sep=""),
                       sep="")
          cat(msg)
          cat("\n")
          invisible(msg)
      })



## ==========================================================================
## curv2Filter
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="curv2Filter"),
          definition=function(object)
      {
          msg <- paste("2D curvature filter '",
                       object@filterId,"' in dimensions ",
                       paste(object@parameters, collapse=" and "),
                       "\nwith settings:",
                       "\n  bwFac=", object@bwFac, "\n  gridsize=",
                       paste(object@gridsize, collapse=",", sep=""),
                       sep="")
          cat(msg)
          cat("\n")
          invisible(msg)
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
          cat("Ellipsoid gate '", identifier(object),
              "' in dimensions ", sep="")
          cat(paste(object@parameters, sep="", collapse=" and "))
          cat("\n")
      })



## ==========================================================================
## kmeansFilter
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="kmeansFilter"),
          definition=function(object)
      {
          msg <- paste("k-means filter '", object@filterId,
                       "' in dimension ", object@parameters[1],
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
          cat(ifelse(length(object@transformation), "transformed", ""),
              "norm2Filter '", identifier(object),
              "' in dimensions ", sep="")
          cat(paste(object@parameters, sep="", collapse=" and "),
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
          nb <-  nrow(object@boundaries)
          cat("Polygonal gate '", identifier(object) ,"' with ",
              ifelse(all(is.na(object@boundaries)), 0, nb),
              " vertices in dimensions ", sep="")
          cat(paste(object@parameters, sep="", collapse=" and "))
          cat("\n")
      })



## ==========================================================================
## quadGate
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
          signature=signature(object="quadGate"),
          definition=function(object)
      {
          cat("Quadrant gate '", identifier(object),
              "' with dimensions:\n", sep="")
          for(i in seq(along=object@parameters)) {
              cat("  ")
              cat(object@parameters[i])
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
          parms <- parameters(object)
          cat("Rectangular gate '", identifier(object),
              "' with dimensions:\n", sep="")
          for(i in seq_along(parms)){
              cat("  ")
              if(is.character(parms[i]))
                  cat(parms[i])
              else
                  cat("anonymous parameter")
              cat(": (")
              cat(paste(object@min[],object@max[i],sep=","))
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
                       "' returning objects with ", object@size," rows", sep="")
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
## workFlow. Recursively print all views in the workflow tree adding
## indentation
## ---------------------------------------------------------------------------
setMethod("show",
          signature=signature(object="workFlow"),
          definition=function(object){
              cat("A flow cytometry workflow called '", object@name,
                  "'\n", sep="")
              tree <- get(object@tree)
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
