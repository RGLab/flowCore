setReplaceMethod("[[",
                 signature(x="filterSet", i="NULL", value="filter"),
                 function(x, i, j, ..., value)
             {
                 n <- paste(".", value@filterId, sep="")
                 x@env[[n]] <- value
                 x@env[[value@filterId]] <- filterReference(x, n)
                 x
             })

## This fucntion was not defined before, and I can only guess that's what
## it was supposed to do initially. No guarantee, though...
.resolveFilter <- function(f)
    get(f@name, f@env)

setReplaceMethod("[[",
                 signature(x="filterSet", i="ANY", value="filterReference"),
                 function(x, i, j, ..., value)
             {
                 ## Copy this filter in from another filterSet
                 if(value@env != x@env)
                     x[[i]] <- .resolveFilter(value)
                 x
             })

setReplaceMethod("[[",
                 signature(x="filterSet", i="NULL", value="formula"),
                 function(x, i, j, ..., value)
             {
                 if(length(value) == 3) {
                     a <- value[[2]]
                     x[[as.character(a)]] <- value
                 } else
                 x[[""]] <- value
                 x
             })

setReplaceMethod("[[",
                 signature(x="filterSet", i="character", value="formula"),
                 function(x, i, j, ..., value)
             {
                 ## Build a place to put names...
                 e <- new.env(parent=x@env)
                 harvest <- function(f) {
                     if(is(f,"name") || is(f,"symbol") ||
                        is(f,"character")) {
                         f <- as.character(f)
                         ## Create new virtual filters to be replaced
                         ## later (or they will throw an error at some
                         ## point).
                         if(!exists(f,env=x@env)) {
                             e[[f]] <- filterReference(x,f)
                         }
                         filterReference(x,f)
                     } else if(length(f) > 1 &&
                               as.character(f[[1]]) %in% c("!","&","|","%on%",
                                                           "%subset%","%&%")) {
                         for(i in 2:length(f)) {
                             f[[i]] <- harvest(f[[i]])
                         }
                         f
                     } else {
                         f
                     }
                 }
                 
                 value <- harvest(value[[length(value)]])
                 ## Extract the names from the formula to resolve to
                 ## formula references
                 old <- parent.env(x@env)
                 parent.env(x@env) <- parent.frame()
                 ## Evaluate in the environment with those names! This
                 ## would be less work if RObjectTables had an R-level
                 ## interface...
                 y <- eval(value,e) 
                 parent.env(x@env) <- old
                 x[[i]] <- y
                 x
             })

setReplaceMethod("[[",
                 signature(x="filterSet", i="missing", value="filter"),
                 function(x, i, j, ..., value) x[[NULL]] <- value)

setReplaceMethod("[[",
                 signature(x="filterSet", i="character", value="filter"),
                 function(x,i,j,...,value)
             {
                 if(nzchar(i))
                     value@filterId <- i
                 x[[NULL]] <- value
                 x
             })

setMethod("[[",
          signature("filterSet", "character"),
          function(x, i, j, ..., drop=TRUE) {
              if(length(i) > 1) stop("Can only select a single item")
              as(x@env[[i]],"concreteFilter")
          })

## Retrieve the filterReferences
setMethod("[",
          signature("filterSet", "character"),
          function(x, i, j, ..., drop=TRUE)
      {
          if(length(i)==1)
              x@env[[i]]
          else
              sapply(i,function(i) x@env[[i]])
      })

setMethod("show",
          signature("filterSet"),
          function(object)
      {
          cat("A set of filter objects:\n")
          cat(paste(sapply(names(object), function(i)
                           identifier(object@env[[i]])), sep=",",
                    collapse=","))
          cat("\n")
      })

setMethod("names",
          signature("filterSet"),
          function(x) ls(env=x@env))

## Performs a topological sort of the filterSet (if possible)
setMethod("sort",
          signature("filterSet"),
          function(x, decreasing=FALSE, dependencies=FALSE, ...)
      {
          n <- names(x)
          ## Dependency matrix for this filter set
          D <- diag(0,length(n))
          row.names(D) <- n
          colnames(D) <- n
          for(i in n) {
              node <- x[[i]]
              if(is(node,"setOperationFilter")) {
                  f <- node@filters
                  ## You can "mask" an operation by building it directly
                  ## into the formula. This supports Gating-ML semantics
                  id <- identifier(node)
                  n2 <- sapply(f,identifier)
                  for(j in f[n2 %in% n]) {
                      D[identifier(j),id] <- 1 
                  }
              }
          }
          ## Parent Set Size
          P <- colSums(D)
          O <- NULL 
          G <- n[P>0]
          while(length(G) > 0) {
              O2 <- n[P==0]
              l <- length(O2)
              if(l > 0) {
                  P[O2] <- -1
                  P <- P - (if(l==1) D[O2,] else apply(D[O2,],2,sum))
                  O <- c(O,O2)
              }
              G <- n[P>0]
              if(length(G)>0 && length(O2)==0)
                  stop(paste("Circular reference detected in filterSet ",
                             "in one or more of:", paste(G, sep=",")))
          }
          ##Add in the leaves
          O <- c(O,n[P==0])
          if(decreasing)
              O <- rev(O)
          if(dependencies)
              attr(O,"AdjM") <- D
          O
      })

setMethod("filterReference",
          signature("filterSet", "character"),
          function(from, name) new("filterReference", env=from@env, name=name)
          )
