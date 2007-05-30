\name{setOperationFilter-class}
\docType{class}
\alias{setOperationFilter-class}

\title{Class "setOperationFilter"}
\description{}

\section{Slots}{
	 \describe{
    \item{\code{filters}:}{Object of class \code{"list"}}
    \item{\code{filterId}:}{Object of class \code{"character"} referencing the filter applied}
    \item{\code{parameters}:}{Object of class \code{"ANY"} describing
      the parameters used to filter the \code{flowFrame}}
    \item{\code{parentId}:}{Object of class \code{"character"}
      referencing the filter previously applied}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filter}"}, directly.
}
\section{Methods}{
No methods defined with class "setOperationFilter" in the signature.
}

\author{B. Ellis}

\seealso{
  \code{\link[flowCore]{filter}}
}
\examples{
showClass("setOperationFilter")
}
\keyword{classes}
