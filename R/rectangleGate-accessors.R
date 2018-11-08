## ==========================================================================
## Methods for objects of type 'rectangleGate'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================






## ==========================================================================
# Compose two rectangle gates together into a higher dimensional cube.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("*",
          signature(e1="rectangleGate",
                    e2="rectangleGate"),
          function(e1, e2)
      {
          if(any(parameters(e1) %in% parameters(e2)))
              stop("Rectangle gate parameters overlap.", call.=FALSE)
          new("rectangleGate", parameters=c(parameters(e1), parameters(e2)),
              min=c(e1@min, e2@min), max=c(e1@max, e2@max))
      })


## ==========================================================================
## subsetting by parameter name 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
setMethod("[",
          signature=signature(x="rectangleGate",
                              i="character"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          if(drop)
              warning("Argument 'drop' ignored for subsetting of ",
                      "rectangelGates", call.=FALSE)
          mt <- i %in% parameters(x)
          if(!all(mt))
              stop("The following parameter(s) are not defined in this gate:",
                   paste("  ", i[!mt], collapse="\n"), call.=FALSE)
          x@min <- x@min[i]
          x@max <- x@max[i]
          parameters(x) <- i
          x
      })


#' @export
setMethod("[",
          signature=signature(x="rectangleGate"),
          definition=function(x, i, j, ..., drop=FALSE)
          stop("rectangleGates may only be subset by parameter name",
               call.=FALSE)
          )
