\name{filterResult-class}
\docType{class}
\alias{filterResult-class}
\alias{==,filterResult,flowFrame-method}
\alias{summary,filterResult-method}
\alias{show,filterResult-method}
\alias{[[,filterResult,ANY-method}
\alias{[[,manyFilterResult,ANY-method}
\alias{[[,multipleFilterResult,ANY-method}
\alias{summary,manyFilterResult-method}
\alias{names,manyFilterResult-method}
\alias{length,manyFilterResult-method}


\title{Class "filterResult" }
\description{Container to store the result of applying a \code{filter} on a
  \code{flowFrame} object}

\section{Slots}{
	 \describe{
    \item{\code{frameId}:}{Object of class \code{"character"}
      referencing the \code{flowFrame} object filtered. Used for sanity checking.}
    \item{\code{filterDetails}:}{Object of class \code{"list"}
      describing the filter applied}
    \item{\code{filterId}:}{Object of class \code{"character"}
      referencing the filter applied}
    \item{\code{parameters}:}{Object of class \code{"ANY"} describing
      the parameters used to filter the \code{flowFrame}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filter}"}, directly.
}
\section{Methods}{
  \item{==}{test equality}
}

\author{B. Ellis, N. LeMeur}

\seealso{
  \code{\link[flowCore]{filter}},
  \code{"\linkS4class{logicalFilterResult}"},
  \code{"\linkS4class{multipleFilterResult}"},
  \code{"\linkS4class{randomFilterResult}"}
}
\examples{
showClass("filterResult")
}
\keyword{classes}