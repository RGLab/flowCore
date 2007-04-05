\name{sampleFilter-class}
\docType{class}
\alias{sampleFilter-class}
\alias{show,sampleFilter-method}

\title{Class "sampleFilter"}
\description{Objects ...}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sampleFilter", ...)}.
	 }
\section{Slots}{
  \describe{
    \item{\code{size}:}{Object of class \code{"numeric"}}
    \item{\code{filterId}:}{A character vector that identifies this
      \code{filter}. This is typically user specified but can be
      automatically deduced by certain filter operations, particularly
      boolean and set operations.}
    \item{\code{parameters}:}{A character vector that identifies the name of
      the parameters that this filter operates upon.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filter}"}, directly.
}
\section{Methods}{
  \item{show}{add description here}
}

\author{B. Ellis}

\seealso{ \code{\linkS4class{filter}} 
}
\examples{
showClass("sampleFilter")
}
\keyword{classes}
