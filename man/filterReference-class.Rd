\name{filterReference-class}
\docType{class}
\alias{filterReference-class}
\alias{filterReference,environment,character-method}


\title{Class filterReference}
\description{A reference to another filter inside a reference. Users should generally
not be aware that they are using this class, but it is used heavily by \code{\linkS4class{filterSet}}
classes.}
\section{Objects from the Class}{
Objects are generally not created by users so there is no constructor function.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{The R name of the referenced filter}
    \item{\code{env}:}{The environment where the filter must live}
    \item{\code{filterId}:}{The filterId, not really used since you always resolve}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filter}"}, directly.
}
\author{ B. Ellis }
\seealso{
	\code{\linkS4class{filterSet}}
}
\examples{
}
\keyword{classes}
