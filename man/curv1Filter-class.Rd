\name{curv1Filter-class}
\docType{class}
\alias{curv1Filter-class}
\alias{curv1Filter} 
\alias{curv1Filter-class}
\alias{show,curv1Filter-method}

\title{Class "curv1Filter" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("curv1Filter", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{bwFac}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{gridsize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"character"} ~~ }
    \item{\code{filterId}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{parameterFilter}"}, directly.
Class \code{"\linkS4class{concreteFilter}"}, by class "parameterFilter", distance 2.
Class \code{"\linkS4class{filter}"}, by class "parameterFilter", distance 3.
}
\section{Methods}{
  \describe{
    \item{\%in\%}{\code{signature(x = "flowFrame", table = "curv1Filter")}: ... }
    \item{show}{\code{signature(object = "curv1Filter")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

 
\seealso{
	}
\examples{
showClass("curv1Filter")
}
\keyword{classes}
