## ==========================================================================
## Methods for objects of type 'ellipsoidGate'
## Note: All filtering methods are stored in file 'in-methods.R'
## ==========================================================================


## ==========================================================================
## show method
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethod("show",
          signature(object="ellipsoidGate"),
          function(object)
      {
          cat("Ellipsoid gate '", identifier(object),
              "' in dimensions ", sep="")
          cat(paste(object@parameters, sep="", collapse=" and "))
          cat("\n")
      })



plotEll <- function(cov, mean, ...)
{
    chol.cov <- t(chol(cov))
    radius <- 1
    t <- seq(0, 2 * base::pi, length = 50)
    ans <- mean +
        (chol.cov %*% rbind(x = radius * cos(t),
                            y = radius * sin(t)))
    lines(t(ans), ...)
}
