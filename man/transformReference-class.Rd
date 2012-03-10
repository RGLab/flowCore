\name{transformReference-class}
\docType{class}
\alias{transformReference-class}
\alias{transformReference}
\alias{parameters,transformReference-method}
\alias{eval,transformReference,missing-method}

\title{Class "transformReference" }
\description{ Class allowing to reference transforms, for instance as parameters.}
\section{Objects from the Class}{
  Objects will be created internally whenever necessary and this should
  not be of any concern to the user.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{The list of references }
    \item{\code{searchEnv}:}{The environment into which the reference points. }
    \item{\code{transformationId}:}{The name of the transformation }
  }
}
\section{Extends}{
Class \code{"\linkS4class{transform}"}, directly.
Class \code{"\linkS4class{transformation}"}, by class "transform", distance 2.
Class \code{"\linkS4class{characterOrTransformation}"}, by class "transform", distance 3.
}
\section{Methods}{
No methods defined with class "transformReference" in the signature.
}

\author{N. Gopalakrishnan}

\keyword{classes}
