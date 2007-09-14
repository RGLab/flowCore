\name{expressionFilter-class}
\docType{class}
\alias{expressionFilter-class}
\alias{expressionFilter}
\alias{show,expressionFilter-method}


\title{Class "expressionFilter"}
\description{A class which holds an arbitrary expression to be used as a filter.}
\section{Objects from the Class}{
Objects can be created with the \code{\link{expressionFilter}} constructor.
}
\section{Slots}{
	 \describe{
    \item{\code{expr}:}{The expression that will be evaluated in the context of the expression values.}
    \item{\code{args}:}{An environment providing additional parameters.}
    \item{\code{filterId}:}{The identifier of the filter}
  }
}
\section{Extends}{
Class \code{"\linkS4class{concreteFilter}"}, directly.
Class \code{"\linkS4class{filter}"}, by class "concreteFilter", distance 2.
}
\author{ B. Ellis }
\seealso{
	\code{\link{expressionFilter}}
}
\keyword{classes}
