\name{lnTransform}
\alias{lnTransform}

\title{Create the definition of a ln transformation function (natural logarthim) to be applied on a data set }
\description{Create the definition of the ln Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x<-log(x)*(r/d).
  The transformation would normally be used to convert to a linear valued parameter
  to the natural logarithm scale.
  Typically r and d are both equal to 1.0. Both must be positive.
}
\usage{
lnTransform(transformationId="defaultLnTransform", r=1, d=1)
}

\arguments{
  \item{transformationId}{character string to identify the transformation}
  \item{r}{positive double that corresponds to a scale factor. }
  \item{d}{positive double that corresponds to a scale factor}
}

\value{
  Returns an object of class \code{transform}.
}

\author{B. Ellis and N. LeMeur}

\seealso{ \code{\link{transform-class}}, \code{\link{transform}}}
\examples{
  data(GvHD)
  lnTrans <- lnTransform(transformationId="ln-transformation", r=1, d=1)
  ln1 <- transform(GvHD,`FSC-H`=lnTrans(`FSC-H`))

opar = par(mfcol=c(2, 1))
plot(density(exprs(GvHD[[1]])[ ,1]), main="Original")
plot(density(exprs(ln1[[1]])[ ,1]), main="Ln Transform")

}
\keyword{methods}