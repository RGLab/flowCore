\name{curv1Filter-class}
\docType{class}
\alias{curv1Filter-class}
\alias{curv1Filter} 
\alias{curv1Filter-class}
\alias{show,curv1Filter-method}

\title{Class "curv1Filter" }

\description{Define a filter that selects high-density regions in one
dimension}

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("curv1Filter", ...)} or using the constructor
\code{curv1Filter}.  }

\section{Slots}{ 
  \describe{

    \item{\code{bwFac}:}{Object of class \code{"numeric"}. The bandwidth
      factor used for smoothing of the density estimate}
    
    \item{\code{gridsize}:}{Object of class \code{"numeric"}. The size
      of the bins used for density estimation }
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameters used to filter the \code{flowFrame} }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}

\section{Extends}{
Class \code{"\linkS4class{parameterFilter}"}, directly.
Class \code{"\linkS4class{concreteFilter}"}, by class "parameterFilter",
distance 2.
Class \code{"\linkS4class{filter}"}, by class "parameterFilter",
distance 3.
}

\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"curv1Filter")}: ... }
    
    \item{show}{\code{signature(object = "curv1Filter")}: ... }
    
  }
}


\author{ Florian Hahne }


\seealso{
  \code{\link[flowCore]{filter}}
  \code{\link[flowCore]{curv2Filter}}
}

\examples{
showClass("curv1Filter")
}

\keyword{classes}
