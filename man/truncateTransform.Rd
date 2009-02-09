\name{truncateTransform}
\alias{truncateTransform}

\title{Create the definition of a trucate tranformation function to be applied on a data set}
\description{Create the definition of the truncate Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x[x<a] <- a.
  Hence, all values less than a are replaced by a. The typical use would be to replace all 
  values less than 1 by 1.
}
\usage{
truncateTransform(transformationId="defaultTruncateTransform", a=1)
}

\arguments{
  \item{transformationId}{character string to identify the transformation}
  \item{a}{double that correponds to the value at which to truncate}
}

\value{
  Returns an object of class \code{transform}.
}

\author{P. Haaland}

\seealso{ \code{\link{transform-class}}, \code{\link{transform}}}
\examples{
samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  truncateTrans <- truncateTransform(transformationId="Truncate-transformation", a=5)
  dataTransform <- transform(samp,`FSC-H`=truncateTrans(`FSC-H`))

}

\keyword{methods }

