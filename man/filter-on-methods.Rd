\name{filter-on-methods}
\docType{methods}
\alias{\%on\%-methods}
\alias{\%on\%}
\alias{\%on\%,filter,transformList-method}
\alias{\%on\%,transformList,flowFrame-method}
\alias{\%on\%,transformList,flowSet-method}
\title{Working with transformed data}

\description{

flowCore provides two methods of working with transformed channels, both of which employ \code{\link{transform}} to construct
the appropriate transform.


}


\title{Methods for Function \%on\% in Package `flowCore'}
\description{
This operator is used to construct a \code{transformFilter} that first applies a \code{transformList} to the data before
applying the \code{filter} operation. You may also apply the operator to a \code{flowFrame} or \code{flowSet} to obtain
transformed values specified in the list.
}
\author{B. Ellis}
\examples{
samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
plot(transform("FSC-H"=log) \%on\% samp)
}
\keyword{methods}

