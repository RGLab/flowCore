\name{logicalFilterResult-class}
\docType{class}
\alias{logicalFilterResult-class}
\alias{logicalFilterResult}
\alias{summary,logicalFilterResult-method}
\alias{names,logicalFilterResult-method}
\alias{length,logicalFilterResult-method}
\alias{[[,logicalFilterResult,ANY-method}

\title{Class "logicalFilterResult"}
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
	\code{\link[flowCore]{filter}}
}
\examples{
showClass("logicalFilterResult")
}
\keyword{classes}