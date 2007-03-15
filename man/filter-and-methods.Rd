\name{filter-and-methods}
\docType{methods}
\alias{\%\&\%-methods}
\alias{\%\&\%}
\alias{\%\&\%,ANY-method}
\alias{\%subset\%,ANY-method}
\alias{\%subset\%}

\title{Take the intersection of two filters with subset semantics}
\description{

The \code{\%subset\%} operator and its alias \code{\%\&\%} are used to construct an intersection filter
with subset semantics. For deterministic filters, such as rectangle, ellipse, and polygonal gates, the 
result will be exactly equivalent to the usual intersection obtained with the \code{\&} boolean operator.

However, the left hand filter is actually evaluated using the subset of events specified by the right hand filter. This
will affect data-driven filters such as \code{\link{norm2Filter}} as well as the calculation of summary statistics returned
by \code{\link{summary}} when applied to \code{\link{filterResult}}
}

\author{B. Ellis}
\keyword{methods}
