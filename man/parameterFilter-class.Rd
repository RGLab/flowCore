\name{parameterFilter-class}
\docType{class}
\alias{parameterFilter-class}
\alias{initialize,parameterFilter-method}


\title{Class "parameterFilter"}
\description{	 A concrete filter that acts on a set of parameters. }
\section{Objects from the Class}{
	\code{parameterFilter} objects are never created directly. This class
	serves as an inheritance point for filters that depends on particular
	parameters.
}
\section{Slots}{
	 \describe{
    \item{\code{parameters}:}{The names of the parameters employed by this filter}
    \item{\code{filterId}:}{The filter identifier}
  }
}
\section{Extends}{
Class \code{"\linkS4class{concreteFilter}"}, directly.
Class \code{"\linkS4class{filter}"}, by class "concreteFilter", distance 2.
}
\section{Methods}{
No methods defined with class "parameterFilter" in the signature.
}

\author{ B. Ellis }
\keyword{classes}
