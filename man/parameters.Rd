\name{parameters}
\docType{methods}
\alias{parameters}
\alias{parameters,flowFrame,missing-method}
\alias{parameters,filter-method}

\title{Methods to retrieve FCS parameters from a flowFrame}
\description{
  Accessor methods for items in the parameters slot of class the
  AnnotatedDataFrame of a flowFrame object or from a filter definition.
  It lists the \code{parameters} and its  values for a flowFrame or a filter.

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

