\name{print.filterSummary}
\docType{methods}
\alias{print.filterSummary}
\alias{summary}
\alias{summary-methods}
\alias{summary,subsetFilter-method}

\title{Summarize the results after applying a filter}

\description{
 Summarize the results after applying a filter on the \code{flowFrame} object.
}
\usage{
summary(object, ...)
}

\arguments{
  \item{object}{Object or list of object of class \code{filterResult},
  \code{logicalFilterResult} or \code{filter}.}
  \item{\dots}{optional arguments for the summary function.}
}

\value{
a print function for a \code{filterSummary} which provides counts and percentages for a 
filtering operation.
}

\author{B. Ellis}

\seealso{\code{filterResult}}
\examples{
data(GvHD)
dead.gate <- norm2Filter("FSC-H","SSC-H",scale.factor=2, filterId="Dead Cells")
Xcell.gate <- norm2Filter("FSC-H","SSC-H",scale.factor=.6, filterId="Selected X Cell")
filter.result <- filter(GvHD[1:4], (Xcell.gate \%&\% !dead.gate))
invisible(sapply(filter.result,function(x) print(summary(x))))
}
\keyword{methods}
\keyword{manip}
