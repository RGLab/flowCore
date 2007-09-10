\name{concreteFilter-class}
\docType{class}
\alias{concreteFilter-class}
\title{Class "concreteFilter"}
\description{The \code{concreteFilter} serves as a base class for all filters that actually 
implement a filtering process. At the moment this includes all filters except 
\code{\linkS4class{filterReference}}, the only non-concrete filter at present.}
\section{Objects from the Class}{
Objects of this class should never be created directly. It serves only as a point of inheritance.
}
\section{Slots}{
	 \describe{
    \item{\code{filterId}:}{The identifier associated with this class. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{filter}"}, directly.
}
\author{ B. Ellis }
\seealso{
	\code{\linkS4class{parameterFilter}}
}
\keyword{classes}
