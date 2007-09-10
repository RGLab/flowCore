\name{parameters}
\docType{methods}
\alias{parameters}
\alias{parameters,flowFrame,missing-method}
\alias{parameters,filter-method}
\alias{parameters,parameterFilter-method}
\alias{parameters,setOperationFilter-method}
\alias{parameters,filterReference-method}
\alias{parameters,parameterTransform-method}

\title{Obtain information about parameters for flow cytometry objects.}
\description{
Many different objects in \code{flowCore} are associated with one or more
parameters. This includes \code{\link[flowCore]{filter}}, \code{\link[flowCore]{flowFrame}}
and \code{\link[flowCore]{parameterTransform}} objects that all either describe or use parameters.

}
\usage{
parameters(object)
}
\arguments{
  \item{object}{Object of class \code{\link[flowCore]{flowFrame}}.}
  }
\value{
 When applied to a \code{flowFrame} objects, the result will be an \code{AnnotatedDataFrame} describing
the parameters recorded by the cytometer. For other objects it will usually return a list of names used
by the object for its calculations.
}
\author{B. Ellis, N. Le Meur}
\examples{
 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 parameters(samp)
 print(samp@parameters@data)
}
\keyword{methods}

