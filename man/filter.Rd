\name{filter}
\alias{filter}
\alias{filter,flowFrame,filter-method}
\alias{filter,flowFrame,rectangleGate}
\alias{filter,flowFrame,polygonGate}
\alias{filter,flowFrame,norm2Filter}
\alias{filter,flowSet,filter-method}
\alias{filter,flowSet,list-method}
\alias{|,filter,filter-method}
\alias{!,filter-method}
\alias{&,filter,filter-method}
\alias{\%subset\%,filter,filter-method}
\alias{draw,filter,flowFrame-method}
\alias{summary,filter-method}
\alias{show,filter-method}
\alias{length,filter-method}


\title{Filter FCS files}
\description{
This function links filter descriptions to a particular set of flow cytometry data allowing for the lightweight calculation
of summary statistics common to flow cytometry analysis.
}

\usage{
filter(x,filter,...)
}
\arguments{
  \item{x}{Object of class \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.}
  \item{filter}{An object of class \code{\link[flowCore]{filter}} or a list of named filters.}
  \item{...}{optional argument}
}
\details{
	The \code{filter} method conceptually links a filter description, represented by a \code{\link[flowCore]{filter}} object,
	to a particular \code{\link[flowCore]{flowFrame}}. This is accomplished via the \code{\link[flowCore]{filterResult}} object,
	which tracks the linked \code{flowFrame} as well as caching the results of the filtering operation itself, allowing for fast
	calculation of certain summary statistics such as the percentage of events accepted by the filter. This method exists chiefly 
	to allow the calculation of these statistics without the need to first \code{\link{Subset}} a \code{flowFrame}, which can be
	quite large.
}

\value{
	A \code{\link[flowCore]{filterResult}} object or a list of \code{\link[flowCore]{filterResult}} objects if \code{x}
	is a \code{\link[flowCore]{flowSet}} object. Note that \code{filterResult} objects are themselves filters, allowing them
	to be used in filter expressions or \code{Subset} operations.
}  

\seealso{\code{\link[flowCore]{filter-class}},
  \code{\link[flowCore::filterResult-class]{filterResult-class}}
}

\author{N. Le Meur}
\examples{
 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 rectGate <- rectangleGate(filterId="nonDebris","FSC-H"=c(200,Inf))
 summary(filter(samp,rectGate))

}
\keyword{methods}  

