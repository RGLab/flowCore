\name{filter-on-methods}
\docType{methods}
\alias{\%on\%-methods}
\alias{\%on\%}
\alias{\%on\%,filter,transformList-method}
\alias{\%on\%,filter,transform-method}
\alias{\%on\%,filter,parameterTransform-method}
\alias{\%on\%,parameterTransform,flowFrame-method}
\alias{\%on\%,transform,flowFrame-method}
\alias{\%on\%,transformList,flowFrame-method}
\alias{\%on\%,transformList,flowSet-method}
\alias{\%on\%,ANY,flowSet-method}


\title{Methods for Function \%on\% in Package `flowCore'}
\description{
This operator is used to construct a \code{transformFilter} that first applies a \code{transformList} to the data before
applying the \code{filter} operation. You may also apply the operator to a \code{flowFrame} or \code{flowSet} to obtain
transformed values specified in the list.
}
\author{B. Ellis}
\examples{
samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
plot(transform("FSC-H"=log, "SSC-H"=log) \%on\% samp)


rectangle <- rectangleGate(filterId="rectangleGateI","FSC-H"=c(4.5, 5.5))
sampFiltered <- filter(samp, rectangle \%on\% transform("FSC-H"=log, "SSC-H"=log))
res <- Subset(samp, sampFiltered)

plot(transform("FSC-H"=log, "SSC-H"=log) \%on\% res)

}
\keyword{methods}

