\name{parameters}
\docType{methods}
\alias{parameters}
\alias{parameters,flowFrame,missing-method}

\title{Methods to retrieve parameters of a flowFrame}
\description{
  Accessor methods for items in the parameters slot (of class
  AnnotatedDataFrame).
  It lists the \code{parameters} and its  values for a flowFrame.

}
\usage{
parameters(object)
}
\arguments{
  \item{object}{Object of class \code{\link[flowCore]{flowFrame}}.}
  }
\author{N. Le Meur}
\examples{
 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 parameters(samp)
 print(samp@parameters@data)
}
\keyword{methods}

