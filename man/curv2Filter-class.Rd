\name{curv2Filter-class}
\docType{class}
\alias{curv2Filter-class}
\alias{curv2Filter} 
\alias{show,curv2Filter-method}

\title{Class "curv2Filter"}

\description{

  A filter that selects high-density regions in two dimension.

}

\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("curv2Filter",
    ...)} or using the constructor \code{curv2Filter}. The constructor
  is the recommended way:

}

\usage{

curv2Filter(x, y, filterId="curv2Filter", bwFac=1.2, gridsize=rep(151, 2))

}

\arguments{

  \item{x,y}{ Characters giving the names of the measurement parameter
    on which the filter is supposed to work on. \code{y} can be missing
    in which case \code{x} is expected to be a character vector of
    length 2 or a list of characters. }
    
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    slot of this gate. }
  
  \item{bwFac, gridsize}{ Numerics of length 1 and 2, respectively, used
    to set the \code{bwFac} and \code{gridsize} slots of the object. }
  
}


\section{Slots}{ 
  \describe{

    \item{\code{bwFac}:}{Object of class \code{"numeric"}. The bandwidth
      factor used for smoothing of the density estimate. }
    
    \item{\code{gridsize}:}{Object of class \code{"numeric"}. The size
      of the bins used for density estimation. }
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameters used to filter the \code{flowFrame}. }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}


\section{Extends}{
  
  Class \code{"\linkS4class{parameterFilter}"}, directly.  Class
  \code{"\linkS4class{concreteFilter}"}, by class
  \code{parameterFilter}, distance 2.  Class
  \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
  distance 3.

}

\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"curv2Filter")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user. }
    
    \item{show}{\code{signature(object = "curv2Filter")}: Print
      information about the filter. }
    
  }
}

\details{

  Areas of high local density in two dimensions are identified by
  detecting significant curvature regions. See \cite{Duong, T. and
  Cowling, A. and Koch, I. and Wand, M.P., Computational Statistics and
  Data Analysis 52/9, 2008} for details.

}


\value{

  Returns a \code{\link{curv2Filter}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}

\author{ Florian Hahne }

\seealso{
  \code{\link[flowCore]{filter}}
  \code{\link[flowCore]{curv1Filter}}
}

\examples{

curv2Filter("FSC-H", "SSC-H", filterId="test filter")

}

\keyword{classes}
\keyword{methods}