\name{singleParameterTransform-class}
\docType{class}
\alias{singleParameterTransform-class}
\alias{initialize,singleParameterTransform-method}

\title{Class "singleParameterTransform"  }

\description{ A transformation that operates on a single parameter}

\section{Objects from the Class}{
  
  Objects can be created by calls of the form
  \code{new("singleParameterTransform", ...)}.
 
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"}. The transformation. }
    \item{\code{parameters}:}{Object of class
      \code{"transformation"}. The parameter to transform. Can be a
      derived parameter from another transformation. }
    \item{\code{transformationId}:}{Object of class
      \code{"character"}. An identifier for the object. }
  }
}

\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}

\section{Methods}{
No methods defined with class "singleParameterTransform" in the signature.
}

\author{F Hahne }

\examples{
showClass("singleParameterTransform")
}

\keyword{classes}
