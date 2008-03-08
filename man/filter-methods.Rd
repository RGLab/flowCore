\name{filter}
\alias{filter}
\alias{filter,flowFrame-method}
\alias{filter,flowFrame,filter-method}
\alias{filter,flowFrame,rectangleGate}
\alias{filter,flowFrame,polygonGate}
\alias{filter,flowFrame,norm2Filter}
\alias{filter,flowFrame,filterSet-method}
\alias{filter,flowSet,filter-method}
\alias{filter,flowSet,filterSet-method}
\alias{filter,flowSet,list-method}
\alias{draw,filter,flowFrame-method}
\alias{summary,filter-method}
\alias{show,filter-method}
\alias{length,filter-method}
\alias{formula,filter-method}
\alias{character,filter-method}
\alias{name,filter-method}
\alias{call,filter-method}
\alias{identifier<-,filter,character-method}

\title{Filter FCS files}

\description{
  
  This function links filter descriptions to a particular set of flow
  cytometry data allowing for the lightweight calculation of summary
  statistics common to flow cytometry analysis.

}

\usage{
filter(x,filter,...)
}
\arguments{
  
  \item{x}{Object of class \code{\linkS4class{flowFrame}} or
    \code{\linkS4class{flowSet}}.}
  
  \item{filter}{An object of class \code{\linkS4class{filter}} or a
    named list \code{filters}.}
  
  \item{...}{Optional arguments}
  
}

\details{
  
  The \code{filter} method conceptually links a filter description,
  represented by a \code{\linkS4class{filter}} object, to a particular
  \code{\linkS4class{flowFrame}}. This is accomplished via the
  \code{\linkS4class{filterResult}} object, which tracks the linked
  frame as well as caching the results of the filtering operation
  itself, allowing for fast calculation of certain summary statistics
  such as the percentage of events accepted by the \code{filter}. This
  method exists chiefly to allow the calculation of these statistics
  without the need to first \code{\link{Subset}} a
  \code{\linkS4class{flowFrame}}, which can be quite large.
	
}

\value{
  
  A \code{\linkS4class{filterResult}} object or a list of
  \code{\linkS4class{filterResult}} objects if \code{x} is a
  \code{\linkS4class{flowSet}}. Note that
  \code{\linkS4class{filterResult}} objects are themselves filters,
  allowing them to be used in filter expressions or \code{Subset}
  operations.

}  

\seealso{\code{\linkS4class{filter}}, \code{\linkS4class{filterResult}}}
}

\author{F Hahne, B. Ellis, N. Le Meur}

\examples{

 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 rectGate <- rectangleGate(filterId="nonDebris","FSC-H"=c(200,Inf))
 summary(filter(samp,rectGate))

}
\keyword{methods}  

