\name{flowSet}
\alias{flowSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a flowSet }
\description{
	A convenience to create a \code{\link{flowSet-class}} object.
}
\usage{
flowSet(..., phenoData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ Any number of \code{flowFrame} objects, either as a list or directly. }
  \item{phenoData}{ An optional \code{\link[Biobase]{AnnotatedDataFrame}} object to attach to the \code{flowSet} }
}
\details{
	This function simply coerces a list to a \code{flowSet} and then attaches the \code{phenoData} parameter
	if present. 
}
\value{
	Returns a \code{flowSet} object.
}
\author{B. Ellis}
\seealso{ \code{\link{flowSet-class}} }
\keyword{methods}
