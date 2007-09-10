\name{sampleFilter}
\alias{sampleFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A filter that selects a random subset of events.}
\description{
This non-parameter filter selects a number of events from the primary \code{flowFrame}.
}
\usage{
sampleFilter(filterId = "sample", size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filterId}{The name of the filter.}
  \item{size}{The number of events to select.}
}
\details{
	Selects a number of events without replacement from a \code{flowFrame}.
}
\value{
	A filter of the class \code{\linkS4class{sampleFilter}}
}
\author{B. Ellis}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(GvHD)
ss   <- sampleFilter(filterId="smallBit", 100)
pop1 <- filter(GvHD,ss)
}
\keyword{methods}
