\name{kmeansFilter-class}
\docType{class}
\alias{kmeansFilter}
\alias{kmeansFilter-class}
\alias{length,kmeansFilter-method}
\alias{show,kmeansFilter-method}

\title{One dimensional k-means filters}
\description{A filter that performs one-dimensional k-means (Lloyd-Max) clustering on a parameter}
\usage{
kmeansFilter(filterId = "kmeans",\dots)
}

\arguments{
  \item{filterId}{ An optional parameter that sets the \code{filterId} of this gate. }
  \item{\dots}{ This gate takes a parameter name and a list of populations }
  
}
\details{
	The one-dimensional k-means filter is a multiple population filter capable of operating on a 
single parameter. It takes a parameter argument associated
with two or more populations. Populations are considered to be ordered such that the population with the 
smallest mean intensity will be the first population in the list and the population with the highest mean
intensity will be the last population listed.
}
\value{
	Returns a \code{\link{kmeansFilter}} object for use in filtering \code{\link{flowFrame}}s or 
	other flow cytometry objects.
}
\author{ B. Ellis and N. LeMeur}
\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }
\examples{
a.1 = kmeansFilter("a","FSC-H"=c("Pop1","Pop2","Pop3"))
summary(filter(read.FCS(system.file("extdata","0877408774.B08", package="flowCore")),a.1))
}
\keyword{methods}

