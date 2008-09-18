\name{parameterTransform-class}
\docType{class}

\alias{parameterTransform-class}

\title{Class "parameterTransform" }
\description{ Link a transformation to particular flow parameters }

\section{Objects from the Class}{

  Objects are created by using the \code{\%on\%} operator and are
  usually not directly instantiated by the user.
  
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"}, the
	 transformation function. }
    \item{\code{parameters}:}{Object of class \code{"character"} The
	 parameters the transformation is applied to.  }
    \item{\code{transformationId}:}{Object of class
	 \code{"character"}. The identifier for the object. }
  }
}

\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{function}"}, by class "transform", distance 2.
}

\section{Methods}{
  \describe{
    \item{\%on\%}{\code{signature(e1 = "filter", e2 =
  "parameterTransform")}: Apply the transformation. }
    \item{\%on\%}{\code{signature(e1 = "parameterTransform", e2 =
  "flowFrame")}: see above }
    \item{parameters}{\code{signature(object = "parameterTransform")}:
  Accessor to the parameters slot }
	 }
}

\author{ Byron Ellis }

\keyword{classes}
