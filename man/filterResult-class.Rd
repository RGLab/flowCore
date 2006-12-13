\name{filterResult-class}
\docType{class}
\alias{filterResult-class}
\alias{filterResult-accessors}
\alias{filterResult}

\title{'filterResult': a filter class representing a filter resolved for a particular data set.}


\description{The \code{filterResult} class describes a realization
 of an abstract filter against a particular \code{\link{flowFrame}}. Being itself a filter
the filter result can participate all other filter compositions. Additionally, \code{filterResult}
objects support a summary method for obtaining simple count statistics relative to the entire data
set. Combined with the \code{\link{\%subset\%}} operator, these summary statistics can also be 
calculated conditional on \code{filterResult} object.
}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{  new("filterResult",}\cr
\code{    frameId=....,   # character string identifying the \code{\link{flowFrame}} used as the data.}
\code{    subSet= ....,   # logical vector}\cr
\code{    filterDetails=...., # a list of details about the filtering operation}\cr
}}

\section{Slots}{
  \describe{
    \item{\code{subSet}:}{Logical \code{vector} identifying data points in the region of interest.}
	\item{\code{frameId}:}{Character \code{vector} identifying the \code{\link{flowFrame}} for this result. Used to check consistency.}
    \item{\code{filterDetails}}{A list giving details about the filtering process.}
  }
}

\section{Methods}{
  \describe{
	\item{\code{summary}:}{Obtain a \code{\link{filterSummary}} object containing simple count statistics.}
    \item{\code{plot:}}{plot the results}
  }
}


\author{B. Ellis, P.D. Haaland, N. Le Meur}

\keyword{classes}