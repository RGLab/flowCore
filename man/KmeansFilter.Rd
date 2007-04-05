\name{kmeansFilter-class}
\docType{class}
\alias{kmeansFilter}
\alias{kmeansFilter-class}
\alias{length,kmeansFilter-method}
\alias{show,kmeansFilter-method}

\title{A multiple population filter determined by kmeans}
\description{
  A convenience method to facilitate the construction of  kmeans \code{\link{filter}} objects
}
\usage{
kmeansFilter(filterId = "kmeans",\dots)
}

\arguments{
  \item{filterId}{ An optional parameter that sets the \code{filterId} of this gate. }
  \item{\dots}{ You can also directly describe a gate without wrapping it in a list, as described below. }
  
}
\details{
  This function only operate on a single parameter.
}
\value{
	Returns a \code{\link{kmeansFilter}} object for use in filtering \code{\link{flowFrame}}s or 
	other flow cytometry objects.
}
\author{ B. Ellis and N. LeMeur}
\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }
\examples{


}
\keyword{methods}

