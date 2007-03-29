\name{logTransform}
\alias{logTransform}

\title{Create the definition of a log transformation function (base specified by user) to be applied on a data set }
\description{Create the definition of the log Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x<-log(x,logbase)*(r/d).
  The transformation would normally be used to convert to a linear valued parameter
  to the natural logarithm scale.
  Typically r and d are both equal to 1.0. Both must be positive.
  logbase = 10 corresponds to base 10 logarithm.
}
\usage{
logTransform(transformationId, logbase=10, r=1, d=1)
}

\arguments{
  \item{transformationId}{character string to identify the transformation}
  \item{logbase}{positive double that correponds to the base of the logarithm. }
  \item{r}{positive double that correponds to a scale factor. }
  \item{d}{positive double that correponds to a scale factor}
}

\value{
  Returns an object of class \code{transform}.
}

\author{B. Ellis}

\seealso{ \code{\link{transform-class}}, \code{\link{transform}}}
\examples{
samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  logTrans <- logTransform(transformationId="ln-transformation", logbase=10, r=1, d=1)
  dataTransform <- transform(samp,`FSC-H`=logTrans(`FSC-H`))
}

\keyword{ methods }