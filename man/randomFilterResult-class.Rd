\name{randomFilterResult-class}
\docType{class}
\alias{randomFilterResult-class}

\title{Class "randomFilterResult"}
\description{Container to store the result of applying a \code{filter} on a
  \code{flowFrame} object}

\section{Slots}{
	 \describe{
    \item{\code{subSet}:}{Object of class \code{"numeric"}}
    \item{\code{frameId}:}{Object of class \code{"character"}  referencing the \code{flowFrame} object filtered. Used for sanity checking.}
    \item{\code{filterDetails}:}{Object of class \code{"list"} describing the filter applied}
    \item{\code{filterId}:}{Object of class \code{"character"} referencing the filter applied}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filterResult}"}, directly.
Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
}

\author{B. Ellis}


\seealso{
	\code{\link[flowCore:filter-methods]{filter}}
}
\keyword{classes}
