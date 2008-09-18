## ==========================================================================
## workFlow objects are basically evaluation environments and an associated
## tree that binds the objects in the environment in a work flow structure.
## All the usual operations for environments (get, assign, ls, rm) are also
## available for workFlow objects.
## ==========================================================================






## ==========================================================================
## Resolve a reference, i.e., get the symbol 'ID' from the evaluation
## environment in 'workspace'.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("get",
          signature=signature(x="character", pos="workFlow",
          envir="missing", mode="missing", inherits="missing"),
          definition=function(x, pos)
      {
          if(!x %in% ls(pos@env)){
              mess <- paste("Unable to resolve reference to object '",
                            x, "'", sep="")
              if(!is.na(w <- pmatch(tolower(x), 
                                    tolower(ls(pos@env))))) 
                  mess <- paste(mess, sprintf("Perhaps you meant %s?",
                            sQuote(ls(pos@env)[w])), sep="\n")
              stop(mess, call.=FALSE)
          }
          get(x, envir=pos@env)
      })

## The same behaviour as above, but allow workflow to be the 'envir' argument
setMethod("get",
          signature=signature(x="character", pos="missing",
          envir="workFlow", mode="missing", inherits="missing"),
          definition=function(x, envir) get(x, envir=envir@env))

## get multiple objects
setMethod("mget",
          signature=signature(x="character", envir="workFlow",
          mode="missing", ifnotfound="missing", inherits="missing"),
          definition=function(x, envir) sapply(x, get, envir))



## ==========================================================================
## List the content of 'env' in a workFlow object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ls",
          signature=signature(name="workFlow", pos="missing",
          envir="missing", all.names="missing", pattern="missing"),
          definition=function(name) ls(name@env))

## Allow for pattern search
setMethod("ls",
          signature=signature(name="workFlow", pos = "missing",
          envir="missing", all.names="missing", pattern="character"),
          definition=function(name, pattern) ls(name@env, pattern=pattern))



## ==========================================================================
## Remove an object from the workFlow environment
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The general method to access the object by name. Note that there there
## are more sophisticated Rm methods for views and actionItems that traverse
## the tree and also remove all dependent objects. Those will get called
## automatically.
setMethod("Rm",
          signature=signature(symbol="character",
                              envir="workFlow",
                              subSymbol="character"),
          definition=function(symbol, envir, subSymbol, ...)
      {
          Rm(get(symbol, envir), envir)
      })



## ==========================================================================
## The names of the views in the workFlow object. Note that this method also
## affects completion for workFlow objects, as only the view are available
## for that. 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("names",
          signature=signature(x="workFlow"),
          definition=function(x){
              nam <- nodes(get(x@tree))
              names(nam) <- sapply(nam, function(y, wf)
                                   identifier(action(get(y, x))), x)
              return(nam)
              
          })


## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## by name
setMethod("[[","workFlow",
          function(x, i, j, ...)
      {
          if(length(i) != 1)
              stop("subscript out of bounds (index must have length 1)",
                   call.=FALSE)
          if(!is.character(i))
              stop("'workFlow' object can only be subset by name",
                   call.=FALSE)
          if(!missing(j))
              warning("Invalid dimension is ignored", call.=FALSE)
          if(! i %in% names(x))
              stop("'", i, "' is not a view in this workFlow object",
                   call.=FALSE)
          get(i,x)
      })

## A useful error message
setMethod("[","workFlow",
          function(x, i, j, ..., drop=FALSE)
      {
          i <- i[[1]]
          stop("This is not a valid operation for a 'workFlow' object.\n",
               "Did you intend to do ", substitute(x), "[[\"",
               substitute(i), "\"]] ?", call.=FALSE)
          })

## By name with completion (but only for views)
setMethod("$",
          signature=signature(x="workFlow",
                              name="character"),
          definition=function(x,name) x@env[[name]])
          


## ==========================================================================
## Plot the workflow tree using Rgraphviz
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Find common action items for each tree level
traverseEdge <- function(g, node=nodes(g)[1], result=NULL)
{
    if(length(node)){
        children <- unlist(adj(g, node))
        lc <- length(children)
        if(lc){
            grps <-
                split(children, sapply(edgeData(g, from=node,
                                                attr="actionItem"),
                                       identifier))
            for(j in grps){
                lg <- length(j)
                result <- rbind(result, c(from=node,
                                          to=j[ceiling(lc/2)]))
            }
            for(i in children)
                result <- traverseEdge(g, i, result=result)
        }
    }
    return(result)
}

setMethod("plot",
          signature=signature(x="workFlow", y="missing"),
          definition=function(x, y, ...)
      {
          if(!suppressWarnings(require(Rgraphviz)))
              stop("You need to have Rgraphviz installed for this feature",
                   call.=FALSE)
          tree <- get(x@tree)
          labels <- gsub("_", "_\n", names(x))
          labels <- paste(sapply(mget(names(x), x), names), labels,
                          sep="\n")
          col <- as.integer(factor(sapply(names(x),
                                          function(y) class(get(y, x)))))
          names(labels) <- names(col) <- names(x)
         
          nodeRenderInfo(tree) <- list(shape="rect", fixedSize=FALSE,
                                       label=labels, col=col, lwd=1,
                                       fontsize=13, textCol=col,
                                       fill="lightgray")
          tmp <- traverseEdge(tree)
          nAttrs <- list()
          if(!is.null(tmp)){
              elabels <-  mapply(function(...)
                                 identifier(edgeData(..., at="actionItem")[[1]]),
                                 from=tmp[,1], to=tmp[,2],
                                 MoreArgs=list(self=tree))
              elabels <-  gsub("_", "_\n", elabels)
              names(elabels) <- apply(tmp, 1, paste, collapse="~")
              edgeRenderInfo(tree) <- list(label=elabels, lwd=2, fontsize=10,
                                       textCol="gray", col="gray")
              width <- rep(1.5, length(names(x)))
              height <- rep(0.8, length(names(x)))
              names(width) <- names(height) <- names(x)
              nAttrs <- list(width=width, height=height)
              g <- Rgraqphviz:::layoutGraph(tree, layoutType="dot", nodeAttrs=nAttrs)  
              Rgraqphviz:::renderGraph(g)
              return(invisible(g))
          }else{
              plot(1, 1, type="n", ann=FALSE, axes=FALSE)
              text(1, 1, sprintf("baseview %s", identifier(x[[names(x)[1]]])))
              return(invisible(tree))
          }
          
         
      })


