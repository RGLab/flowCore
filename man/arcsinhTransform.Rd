\name{arcsinhTransform}
\alias{arcsinhTransform}

\title{Create the definition of an arcsinh transformation function (base specified by user) to be applied on a data set }
\description{Create the definition of the arcsinh Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x<-asinh(a+b*x)+c).
  The transformation would normally be used to convert to a linear valued parameter
  to the natural logarithm scale. By default a and b are both equal to 1 and c to 0. 
  
}
\usage{
arcsinhTransform(transformationId="defaultArcsinhTransform", a=1, b=1, c=0)
}

\arguments{
  \item{transformationId}{character string to identify the transformation}
  \item{a}{positive double that corresponds to a shift about 0. }
  \item{b}{positive double that corresponds to a scale factor.}
  \item{c}{positive double}
}

\value{
  Returns an object of class \code{transform}.
}

\author{B. Ellis}

\seealso{ \code{\link{transform-class}},
  \code{\link{transform}}, \code{asinh}
}
\examples{
samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  asinhTrans <- arcsinhTransform(transformationId="ln-transformation", a=1, b=1, c=1)
  translist <- transformList('FSC-H', asinhTrans) 
  dataTransform <- transform(samp, translist)
}

\keyword{ methods }
