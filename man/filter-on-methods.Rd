\name{filter-on-methods}
\docType{methods}
\alias{\%on\%-methods}
\alias{\%on\%}
\alias{\%on\%, filter,transformList}
\alias{\%on\%, transformList,flowFrame}
\alias{\%on\%, transformList, flowSet}


\title{Methods for Function \%on\% in Package `flowCore'}
\description{
This operator is used to construct a \code{transformFilter} that first applies a \code{transformList} to the data before
applying the \code{filter} operation. You may also apply the operator to a \code{flowFrame} or \code{flowSet} to obtain
transformed values specified in the list.
}
\author{B. Ellis}
\examples{
\dontrun{
samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
plot(transform("FSC-H"=log) %on% samp)
}
}
\keyword{methods}

