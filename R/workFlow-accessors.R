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
              envir="missing"),
          definition=function(x, pos=-1, envir=as.environment(pos),
              mode="any", inherits=TRUE)
      {
          checkClass(x, "character", 1)
          allNames <- c(ls(pos), ls(alias(pos)))
          if (!x %in% allNames) {
              mess <- paste("Unable to resolve reference to object '",
                            x, "'", sep = "")
              if (!is.na(w <- pmatch(tolower(x), tolower(allNames))))
                  mess <- paste(mess, sprintf("Perhaps you meant %s?",
                                              sQuote(ls(pos@env)[w])), sep = "\n")
              stop(mess, call. = FALSE)
          }
          id <- if (!x %in% ls(pos))
              alias(pos)[[x]]
          else x
          if (length(id) != 1)
              stop("The alias '", x, "' is not unique.\n Unable to resolve",
                   " to ID.", call. = FALSE)
          get(id, envir = pos@env)
      })

## The same behaviour as above, but allow workflow to be the 'envir' argument
setMethod("get",
          signature=signature(x="character", pos="missing",
              envir="workFlow"),
          definition=function(x, pos=-1, envir=as.environment(pos),
              mode="any", inherits=TRUE) get(x, envir=envir@env))

## get multiple objects
setMethod("mget",
          signature=signature(x="character", envir="workFlow"),
          definition=function(x, envir, mode="any",
              ifnotfound=list(function(x) stop(paste0("value for '",
                  x, "'not found"), call.=FALSE))) sapply(x, get, envir))



## ==========================================================================
## List the content of 'env' in a workFlow object
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("ls",
          signature=signature(name="workFlow",
                              pos="missing",
          envir="missing", all.names="missing", pattern="missing"),
          definition=function(name) ls(name@env))

## Allow for pattern search
setMethod("ls",
          signature=signature(name="workFlow",
                              pos = "missing",
          envir="missing", all.names="missing", pattern="character"),
          definition=function(name, pattern) ls(name@env, pattern=pattern))



## ==========================================================================
## Subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## is a view name unique?
isUniqueName <- function(n, wf, warn=FALSE, uid)
{
    atab <-  mget(ls(alias(wf)), alias(wf))
    if(length(atab[[n]])==0)
            stop("'", n, "' is not a view in this workFlow object.",
                 call.=FALSE)
    res <- length(atab[[n]])==1
    if(!res && warn)
        warning("'", n, "' is not a unique view name in this workFlow object.\n",
                "You will only be able to access the view using its UID",
                " '", uid, "'", call.=FALSE)
     return(res)
}

## by name
setMethod("[[",
          signature=signature(x="workFlow"),
          definition=function(x, i, j, ...)
      {
          if(length(i) != 1)
              stop("subscript out of bounds (index must have length 1)",
                   call.=FALSE)
          if(!is.character(i))
              stop("'workFlow' object can only be subset by name",
                   call.=FALSE)
          if(!missing(j))
              warning("Invalid dimension is ignored", call.=FALSE)
          ids <- nodes(tree(x))
          if(i %in% ids)
              return(get(i, x))
          if(!isUniqueName(i, x))
              stop("'", i, "' is not a unique name for a view in this workflow.\nIn",
                   " order to access it you have to use its unique identifier. ",
                   "\n'", i, "' maps to the following UIDs:\n  '",
                   paste(alias(x)[[i]], collapse="', '"), "'", call.=FALSE)
          get(i,x)
      })

## A useful error message
setMethod("[",
          signature=signature(x="workFlow"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          i <- i[[1]]
          stop("This is not a valid operation for a 'workFlow' object.\n",
               "Did you intend to do ", substitute(x), "[[\"",
               substitute(i), "\"]] ?", call.=FALSE)
          })

## By name with completion (but only for views)
setMethod("$",
          signature=signature(x="workFlow"),
          definition=function(x,name) x[[name]])



## ==========================================================================
## Get the alias table from a workFlow or an environment
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
uniqueAlias <- function(x, wf) length(alias(wf)[[x]])==1

setMethod("alias",
          signature=signature(object="workFlow"),
          definition=function(object) get(object@alias))

setMethod("alias",
          signature=signature(object="environment"),
          definition=function(object){
              alRef <- grep("aliasRef", ls(object), value=TRUE)
              if(length(alRef) != 1)
                  stop("Unable to find alias table in this environment.",
                       call.=FALSE)
              get(alRef, object)
          })



## ==========================================================================
## Get the journal from a workFlow or an environment
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("journal",
          signature=signature(object="workFlow"),
          definition=function(object) get(object@journal))

setMethod("journal",
          signature=signature(object="environment"),
          definition=function(object){
              jRef <- grep("journalRef", ls(object), value=TRUE)
              if(length(jRef) != 1)
                  stop("Unable to find journal in this environment.",
                       call.=FALSE)
              get(jRef, object)
          })



## ==========================================================================
## Get the nodes from a workFlow tree
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("nodes",
          signature=signature(object="workFlow"),
          definition=function (object, ...)
      {
          n <- nodes(tree(object))
          names(n) <- id2Alias(n, object)
          n
      })


nameIdMapping <- function(wf)
{
    n <- nodes(wf)
    vec <- names(n)
    names(vec) <- n
    vec
}



## ==========================================================================
## Get the tree from a workFlow or an environment
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("tree",
          signature=signature(object="workFlow"),
          definition=function(object) get(object@tree))

setMethod("tree",
          signature=signature(object="environment"),
          definition=function(object){
              tRef <- grep("treeRef", ls(object), value=TRUE)
              if(length(tRef) != 1)
                  stop("Unable to find tree in this environment.",
                       call.=FALSE)
              get(tRef, object)
          })




## ==========================================================================
## Plot the workflow tree using Rgraphviz
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("plot",
          signature=signature(x="workFlow",
                              y="missing"),
          definition=function(x, y, ...)
      {
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
                      if(length(grps)>2)
                          subGraphs[["node"]] <<- grps
                      for(j in grps){
                          lg <- length(j)
                          result <- rbind(result, c(from=node,
                                                    to=as.vector(j[ceiling(lg/2)])))
                      }
                      for(i in children)
                          result <- traverseEdge(g, i, result=result)
                  }
              }
              return(result)
          }
          if (!suppressWarnings(requireNamespace(Rgraphviz)))
              stop("You need to have Rgraphviz installed for this feature",
                   call.=FALSE)
          tree <- get(x@tree)
          labels <- gsub(" ", "\n", gsub("_", "_\n", views(x)))
          cmatch <- cbind(c("view", "compensateView", "transformView",
                            "gateView"), c("black", "green", "blue", "red"))
          mt <- match(sapply(views(x), function(y) class(get(y, x))),
                      cmatch)
          col <- cmatch[mt, 2]
          nn <- views(x)
          names(labels) <- names(col) <- getAlias(nn, x)
          nodeRenderInfo(tree) <- list(shape="rect", fixedSize=FALSE,
                                       label=labels, col=col, lwd=1, fontsize=16,
                                       textCol=col, fill="lightgray")
          subGraphs <- list()
          tmp <- traverseEdge(tree)
          nAttrs <- list()
          if (!is.null(tmp)) {
              elabels <- id2Alias(mapply(function(...)
                                         identifier(edgeData(..., attr="actionItem")[[1]]),
                                         from=tmp[, 1], to=tmp[, 2],
                                         MoreArgs=list(self=tree)), x)
              elabels <- gsub("_", "_\n", elabels)
              names(elabels) <- apply(tmp, 1, paste, collapse="~")
              edgeRenderInfo(tree) <- list(label=elabels, lwd=2,
                                           fontsize=10, textCol="gray", col="gray")
              width <- rep(1.5, length(nn))
              height <- rep(0.8, length(nn))
              names(width) <- names(height) <- nn
              nAttrs <- list(width=width, height=height)
              g <- Rgraphviz::layoutGraph(tree, layoutType="dot",
                                           nodeAttrs=nAttrs)
              Rgraphviz::renderGraph(g)
              return(invisible(g))
          }
          else {
              plot(1, 1, type="n", ann=FALSE, axes=FALSE)
              text(1, 1, sprintf("baseview %s", identifier(x[[nn[1]]])))
              return(invisible(tree))
          }
      })



## nodePositions <- function(wf)
## {
##     ## some global variables to be set by the recursive function
##     tree <- get(wf@tree)
##     yoffset <- 1
##     xoffset <- 1
##     width <- 0
##     depth <- 0
##     nodes <- matrix(ncol=5, nrow=length(nodes(tree)),
##                     dimnames=list(nodes(tree),
##                                   c("x", "y", "width", "shift", "mids")))
##     ## return all children of a node
##     children <- function(node, tree) unlist(adj(tree, node))
##     ## Called a node is visited the first time in an Euler Tour
##     firstVisit <- function(node)
##     {
##         nodes[node, "x"] <<- width
##         nodes[node, "y"] <<- depth
##         depth <<- depth+yoffset
##     }
##     ## Called a node is visited the last time in an Euler Tour
##     lastVisit <- function(node)
##     {
##         ##textWidth <- max(strwidth(id2Alias(node, wf)), xoffset)
##         textWidth <- xoffset
##         shift <- 0
##         x <- nodes[node, "x"]
##         boxWidth <- width-x
##         if(textWidth > boxWidth){
##             delta <- textWidth-boxWidth
##             boxWidth <- textWidth
##             width <<- width+delta
##             shift <- delta/2
##         }
##         nodes[node, "width"] <<- boxWidth
##         nodes[node, "mids"] <<- as.integer(x+(boxWidth/2))
##         nodes[node, "shift"] <<- shift
##         depth <<- depth-yoffset
##     }
##     ## Called if the node is a leaf node
##     externalVisit <- function(node)
##     {
##         ## textWidth <- max(strwidth(id2Alias(node, wf)), xoffset)
##         textWidth <- xoffset
##         nodes[node, "x"] <<- width
##         nodes[node, "y"] <<- depth
##         nodes[node, "width"] <<- textWidth
##         nodes[node, "mids"] <<- as.integer(width+(textWidth/2))
##         width <<- width+textWidth
##         nodes[node, "shift"] <<- 0
##     }
##     ## recursive function to calculate the bounding box for each node
##     boundBox <- function(g, node=nodes(g)[1])
##     {
##         child <- children(node, g)
##         lc <- length(child)
##         if(!lc) externalVisit(node)
##         firstVisit(node)
##         for(i in child)
##             boundBox(g, i)
##         lastVisit(node)
##     }
##     boundBox(tree)
##     nodes[, "x"] <- nodes[, "x"] - nodes[, "shift"]
##     nodes[, "shift"] <- nodes[, "x"] + nodes[, "width"]
##     nodes[, "y"] <- abs(nodes[, "y"] - max(nodes[, "y"]))
##     colnames(nodes)[c(1,4)] <- c("x1", "x2")
##     rownames(nodes) <- id2Alias(rownames(nodes), wf)
##     return(nodes)
## }


## nodes <- nodePositions(wf)
## plot(1,1, xlim=c(0, max(apply(nodes[, c(1,3)], 1, sum))),
##      ylim=c(0, max(nodes[,"y"])+1), type="n", axes=FALSE, ann=FALSE)
## #plot(1,1, xlim=0:1, ylim=0:1, type="n", axes=FALSE, ann=FALSE)
## par(mar=rep(0.5,4))
## #rect(nodes[,"x1"], nodes[,"y"],
## #     nodes[,"x2"], nodes[,"y"]+1,
## #     col=seq_len(nrow(nodes))+1)
## #rect(nodes[,"x"], nodes[,"y"],
## #     nodes[,"x"]+nodes[,"width"], nodes[,"y"]+1,
## #     col=seq_len(nrow(nodes)))
## cex <- 0.5
## xmids <- nodes[,"x1"]+(nodes[,"width"])/2
## ymids <- nodes[,"y"]+0.5
## #text(xmids, nodes[,"y"]+0.5, rownames(nodes), cex=cex)
## sw <- strwidth(rownames(nodes), cex=cex)*1.2
## sh <- strheight(rownames(nodes), cex=cex)*5
## #rect(xmids-(sw/2), ymids+(sh/2),
## #     xmids+(sw/2), ymids-(sh/2))
## points(xmids, ymids, col="#00000090", pch=20, cex=2)
## abline(v=0:100, col="#00000050")
