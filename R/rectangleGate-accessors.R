## ==========================================================================
## Methods for objects of type 'rectangleGate'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",signature(object="rectangleGate"),function(object)
      {
          parms <- parameters(object)
          cat("Rectangular gate '", identifier(object),
              "' with dimensions:\n", sep="")
          for(i seq_along(parms)) {
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
# Compose two rectangle gates together into a higher dimensional cube.
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("*",
          signature(e1="rectangleGate", e2="rectangleGate"),
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
setMethod("[",
          signature=signature("rectangleGate", "character"),
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
          x@parameters <- i
          x
      })


setMethod("[",
          signature=signature("rectangleGate"),
          definition=function(x, i, j, ..., drop=FALSE)
          stop("rectangleGates may only be subset by parameter name",
               call.=FALSE)
          )
