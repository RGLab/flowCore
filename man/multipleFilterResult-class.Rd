\name{multipleFilterResult-class}
\docType{class}
\alias{multipleFilterResult-class}
\alias{multipleFilterResult}
\alias{length,multipleFilterResult-method}
\alias{names,multipleFilterResult-method}
\alias{names<-,multipleFilterResult-method}
\alias{[[,multipleFilterResult-method}
\alias{[[,multipleFilterResult,ANY-method}
\alias{[,multipleFilterResult,ANY-method}
\alias{summary,multipleFilterResult-method}
\alias{show,multipleFilterResult-method}

\title{Class "multipleFilterResult"}
\description{Container to store the result of applying \code{filter} on
  set of \code{flowFrame} objects}

\section{Slots}{
	 \describe{
    \item{\code{subSet}:}{Object of class \code{"factor"}}
    \item{\code{frameId}:}{Object of class \code{"character"} referencing the \code{flowFrame} object filtered. Used for sanity checking.}
    \item{\code{filterDetails}:}{Object of class \code{"list"} describing the filter applied}
    \item{\code{filterId}:}{Object of class \code{"character"} referencing the filter applied}
    \item{\code{parameters}:}{Object of class \code{"ANY"} describing
      the parameters used to filter the \code{flowFrame}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filterResult}"}, directly.
Class \code{"\linkS4class{filter}"}, by class "filterResult", distance 2.
}
\section{Methods}{
  \item{[, [[}{subsetting. If \code{x} is \code{multipleFilterResult},
    then \code{x[[i]]} a \code{FilterResult}  object. The semantics is
    similar to the behavior of the subsetting operators for lists.}
  \item{length}{number of \code{FilterResult} objects in the set.}
  \item{names}{names of the  \code{FilterResult} objects in the set.}
  \item{summary}{summary \code{FilterResult} objects in the set.}
}

\author{B. Ellis}

\seealso{
	 \code{\link[flowCore]{FilterResult}}
	
}
\examples{
showClass("multipleFilterResult")
}
\keyword{classes}
