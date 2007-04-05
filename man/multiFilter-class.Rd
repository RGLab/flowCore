\name{multiFilter-class}
\docType{class}
\alias{multiFilter-class}

\title{Class to represent list of filter}
\description{Class to represent list of filter}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("multiFilter", ...)}.
	  
}
\section{Slots}{
	 \describe{
    \item{\code{populations}:}{A \code{"character"} vector to identify
      the collection of filters.}
    \item{\code{filters}:}{Object of class \code{"list"}. Each element
      of the list corresponds to a filter}
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
No methods defined with class "multiFilter" in the signature.
}

\author{B. Ellis}

\seealso{
  \code{\linkS4class{filter}}
}
\examples{
showClass("multiFilter")

R1 <- rectangleGate(filterId="rectangleGateI","FSC-H"=c(1,50),"SSC-H"=c(50,800))
R2 <- rectangleGate(filterId="rectangleGateII",list("FSC-H"=c(1,50),"SSC-H"=c(50,800)))

fil <- new("multiFilter",populations = "Pop1", filters=c(R1,R2),
filterId = "filterSchemaI", parameters = unique(c(R1@parameters, R1@parameters)))
}
\keyword{classes}
