\name{scaleTransform}
\alias{scaleTransform}

\title{Create the definition of a scale transformation function to be applied on a data set }
\description{Create the definition of the scale Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x = (x-a)/(b-a).
  The transformation would normally be used to convert to a 0-1 scale. In this case,
  b would be the maximum possible value and a would be the minimum possible value.
}
\usage{
scaleTransform(transformationId, a, b)
}

\arguments{
  \item{transformationId}{character string to identify the transformation}
  \item{a}{double that correponds to the value that will be transformed to 0}
  \item{b}{double that correponds to the value that will be transformed to 1}
}

\value{
  Returns an object of class \code{transform}.
}

\author{P. Haaland}

\seealso{ \code{\link{transform-class}}, \code{\link{transform}}}
\examples{
\dontrun{samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  scaleTrans <- scaleTransform(transformationId="Truncate-transformation", a=1, b=10^4)
  dataTransform <- transform(samp,`FLH-1`=scaleTrans(`FLH-1`))
}
}

\keyword{ methods }

