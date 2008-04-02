\name{curv2Filter-class}
\docType{class}
\alias{curv2Filter-class}
\alias{curv2Filter} 
\alias{show,curv2Filter-method}

\title{Class "curv2Filter"}

\description{A filter that selects high-density regions in tow
  dimension}

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("curv2Filter", ...)} or using the constructor
\code{curv2Filter}.  }

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
	"curv2Filter")}: ... }
    
    \item{show}{\code{signature(object = "curv2Filter")}: ... }
    
  }
}

\author{ Florian Hahne }

\seealso{
  \code{\link[flowCore]{filter}}
  \code{\link[flowCore]{curv1Filter}}
}

\examples{
showClass("curv2Filter")
}

\keyword{classes}
