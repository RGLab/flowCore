\name{logicleTransform}
\alias{logicleTransform}

\title{Compute a transform using the 'logicle' function}
\description{
  Logicle transformations create a subset of \code{biexponential}
  transformation
  logicle(x)= 
 
}
\usage{
logicleTransform(transformationId, w = 0, r = 262144, d = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used by the transform/filter integration routines. }
  \item{w}{positive number that corresponds to the width of the negative data range and the range of linearized
    datain natural log}
  \item{r}{top of the scale data value, e.g, 10000 for common 4 decade
    data or 262144 for a 18 bit data range.}
  \item{d}{breath of the display in natural logarithm units.}
  \item{\dots}{additional arguments}
}
\references{Parks D.R., Roederer M., Moore W.A.(2006)  A new "Logicle" display
  method avoids deceptive effects of logarithmic scaling for low signals
and compensated data. CytometryA, 96(6):541-51.}
\author{ B. Ellis}

\seealso{\code{\link[flowCore]{biexponential}} }
\examples{
\dontrun{
# Construct some "flow-like" data which tends to be hetereoscedastic.
before = sapply(10^sample(c(1,2,4),1000),function(i) rnorm(i,i*1.2))

#Just use the standard arguments
logicle  = LogicleTransform(w=10)
after.1 = logicle(before)
}
  }
}
\keyword{methods}
