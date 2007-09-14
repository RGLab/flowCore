\name{filter-and-methods}
\docType{methods}
\alias{intersectFilter-method}
\alias{subsetFilter-method}
\alias{\%&\%}
\alias{\%&\%-methods}
\alias{\%&\%,ANY-method}
\alias{\%subset\%,ANY-method}
\alias{\%subset\%}
\alias{&,filter,filter-method}
\alias{&,filter,list-method}
\alias{&,list,filter-method}
\alias{\%subset\%,filter,filter-method}
\alias{\%subset\%,list,filter-method}
\alias{coerce,intersectFilter,call-method}

\title{Take the intersection of two filters}
\description{

There are two notions of intersection in \code{flowCore}. First, there is the usual intersection boolean operator
\code{&} that has been overridden to allow the intersection of two filters or of a filter and a list for convenience. There
is also the \code{\%&\%} or \code{\%subset\%} operator that takes an intersection, but with subset semantics rather than simple
intersection semantics. In other words, when taking a subset, calculations from \code{\link[flowCore]{summary}} and other methods are taken
with respect to the right hand filter. This primarily affects percentage calculations, which are ordinarily calculated with respect
to the entire population as well as data-driven gating procedures which will operate only on elements contained by the right hand filter.
This becomes especially important when using filters such as \code{\link[flowCore]{norm2Filter}}

}

\author{B. Ellis}
\keyword{methods}
