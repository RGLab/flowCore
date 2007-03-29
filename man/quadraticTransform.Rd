\name{quadraticTransform}
\alias{quadraticTransform}

\title{Create the definition of a quadratic transformation function to be applied on a data set }
\description{Create the definition of the quadratic Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x <- a*x\^2 + b*x + c
}
\usage{
quadraticTransform(transformationId, a = 1, b = 1, c = 0)
}

\arguments{
  \item{transformationId}{character string to identify the transformation }
  \item{a}{double that correponds to the quadratic coefficient in the equation}
  \item{b}{double that correponds to the linear coefficient in the equation}
  \item{c}{double that correponds to the intercept in the equation}
}
\value{
  Returns an object of class \code{transform}.
}

\author{ N. Le Meur }

\seealso{ \code{\link{transform-class}}, \code{\link{transform}}}
\examples{
samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  quadTrans <- quadraticTransform(transformationId="Quadratic-transformation", a=1, b=1, c=0)
  dataTransform <- transform(samp,`FSC-H`=quadTrans(`FSC-H`))

}

\keyword{methods}
