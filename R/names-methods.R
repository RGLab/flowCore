## ==========================================================================
## Show the name of an object, or of the items contained in a (usually
## list-like) object. In flowCore, names are often used for subpopulations.
## Note that parameter names are accessed by the parameters method.
## ==========================================================================






## ==========================================================================
## Accessor to name slot, a human-readable identifier of the object.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="actionItem"),
          definition=function(x) x@name)



## ==========================================================================
## The identifiers of the flowFrames (i.e., the sampleNames of the flowSet)
## for which the filterResults have been computed. 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="filterResultList"),
          definition=function(x) x@frameId)



## ==========================================================================
## The identifiers of the individual filters in the flowSet
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="filterSet"),
          definition=function(x) ls(env=x@env))



## ==========================================================================
## The names of the populations if the summary was computed for a
## multipleFilterResult, for a logicalFilterResult the name of the input
## filter.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="filterSummary"),
          definition=function(x) x@name)



## ==========================================================================
## This return a pretified version of the parameter names, including the
## parameter description if present
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="flowFrame"),
          definition=function(x)
      {
          cn <- colnames(x)
          fn <- featureNames(x)
          if(length(fn) == length(cn)) {
              cn <- paste("<", cn, ">", sep="")
              for(i in seq(along=fn)) {
                  if(!is.na(fn[i]) && fn[i]!="")
                      cn[i] <- paste(cn[i],fn[i])
              }
          }
          cn
      })



## ==========================================================================
## The names of the two populations created by the filter. For
## logicalFilterResults we always have one population in the filter
## and the complement. 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="logicalFilterResult"),
          definition=function(x) paste(x@filterId, c("+", "-"), sep=""))



## ==========================================================================
## The names of the individual sub-populations (i.e., the colnames of the
## subSet matrix.
## ---------------------------------------------------------------------------
setMethod("names",
          signature=signature(x="manyFilterResult"),
          definition=function(x) colnames(x@subSet))



## ==========================================================================
## The names of the individual sub-populations (i.e., the levels of the
## subSet factor). The replacement method simply changes the factor levels.
## ---------------------------------------------------------------------------
setMethod("names",
          signature=signature(x="multipleFilterResult"),
          definition=function(x) levels(x@subSet))

setReplaceMethod("names",
                 signature=signature(x="multipleFilterResult",
                                     value="ANY"),
                 definition=function(x, value)
             {
                 if(length(value) != length(levels(x@subSet)))
                     stop("Length of replacement vector doesn't match.")
                 levels(x@subSet) <- value
                 x@filterDetails[[1]]$populations <- value
                 return(x)
             })



## ==========================================================================
## Accessor to name slot, a human-readable identifier of the view.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="view"),
          definition=function(x) x@name)



## ==========================================================================
## The names of only the views and actionItems in the workFlow object. Note
## that this method also affects completion for workFlow objects, as only
## view and actionItem references are being completed. Use 'views' or
## 'action', respectively to only get on of the type. 'ls' will give you
## all the symbols from the environment.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Find alias for an identifier
id2Alias <- function(id, workflow)
{
    checkClass(id, "character")
    checkClass(workflow, "workFlow")
    workflow <- alias(workflow)
    fun <- function(y){
        ind <- names(which(sapply(as.list(workflow), function(x)
                                  y %in% x)==TRUE))
        if(length(ind)==1 && length(workflow[[ind]])==1)
            ind
        else
            y
    }
    as.vector(sapply(id, fun))
}

setMethod("names",
          signature=signature(x="workFlow"),
          definition=function(x)
      {
          nam <- nodes(get(x@tree))
          acts <- unique(unlist(sapply(nam, function(y)
                                       identifier(action(get(y, x))))))
          return(id2Alias(c(nam, acts), x))
          
      })



## ==========================================================================
## The names of only the views in the workFlow object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("views",
          signature=signature(x="workFlow"),
          definition=function(x)
          return(id2Alias(nodes(get(x@tree)), x)))    



## ==========================================================================
## The names of only the actionItems in the workFlow object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("actions",
          signature=signature(x="workFlow"),
          definition=function(x)
      {
          nam <- nodes(get(x@tree))
          acts <- unique(unlist(sapply(nam, function(y)
                                       identifier(action(get(y, x))))))
          return(id2Alias(acts, x))
      })
