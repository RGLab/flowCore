\name{manyFilterResult-class}
\docType{class}
\alias{manyFilterResult-class}
\alias{length,manyFilterResult-method}
\alias{names,manyFilterResult-method}
\alias{[[,manyFilterResult-method}
\alias{[[,manyFilterResult,ANY-method}
\alias{summary,manyFilterResult-method}
\alias{show,manyFilterResult-method}
\alias{as.data.frame.manyFilterResult}
\alias{manyFilterResult}
\alias{parameters,manyFilterResult-method}

\title{Class "manyFilterResult"}
\description{The result of a several related, but possibly
overlapping filter results. The usual creator of this object will usually
be a \code{\link{filter}} operation of \code{\link{filterSet}} object on a 
\code{\link{flowFrame}} object.} 

\section{Slots}{
	 \describe{
    \item{\code{subSet}:}{Object of class \code{"matrix"}}
    \item{\code{frameId}:}{Object of class \code{"character"} referencing the \code{flowFrame} object filtered. Used for sanity checking.}
    \item{\code{filterDetails}:}{Object of class \code{"list"} describing the filter applied}
    \item{\code{filterId}:}{Object of class \code{"character"} referencing the filter applied}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filterResult}"}, directly.
Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
}
\section{Methods}{
  \item{[, [[}{subsetting. If \code{x} is \code{manyFilterResult},
    then \code{x[[i]]} a \code{FilterResult}  object. The semantics is
    similar to the behavior of the subsetting operators for lists.}
  \item{length}{number of \code{filterResult} objects in the set.}
  \item{names}{names of the  \code{filterResult} objects in the set.}
  \item{summary}{summary \code{filterResult} objects in the set.}
}

\author{B. Ellis}

\seealso{
	 \code{\link[flowCore]{FilterResult}}
	
}
\examples{
showClass("manyFilterResult")
}
\keyword{classes}
