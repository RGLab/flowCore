\name{curv2Filter-class}
\docType{class}
\alias{curv2Filter-class}
\alias{curv2Filter} 
\alias{show,curv2Filter-method}


\title{Class "curv2Filter"}


\description{

  Class and constructor for data-driven \code{\link{filter}} objects
  that selects high-density regions in two dimensions.

}


\section{Extends}{
  
  Class \code{"\linkS4class{parameterFilter}"}, directly.

  Class \code{"\linkS4class{concreteFilter}"}, by class
  \code{parameterFilter}, distance 2.

  Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
  distance 3.

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


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("curv2Filter",
    ...)} or using the constructor \code{curv2Filter}. The constructor
  is the recommended way of object instantiation:

}


\usage{

curv2Filter(x, y, filterId="defaultCurv2Filter", bwFac=1.2,
gridsize=rep(151, 2))

}


\arguments{

  \item{x,y}{ Characters giving the names of the measurement parameter
    on which the filter is supposed to work on. \code{y} can be missing
    in which case \code{x} is expected to be a character vector of
    length 2 or a list of characters. }
    
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    slot of this filter. The object can later be identified by this
    name.}
  
  \item{bwFac, gridsize}{ Numerics of length 1 and 2, respectively, used
    to set the \code{bwFac} and \code{gridsize} slots of the object. }
  
}


\value{

  Returns a \code{\link{curv2Filter}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"curv2Filter")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "curv2Filter")}: Print
      information about the filter. }
    
  }
}


\details{

  Areas of high local density in two dimensions are identified by
  detecting significant curvature regions. See \cite{Duong, T. and
  Cowling, A. and Koch, I. and Wand, M.P., Computational Statistics and
  Data Analysis 52/9, 2008} for details. The constructor
  \code{curv2Filter} is a conveniance function for object
  instantiation. Evaluating a \code{curv2Filter} results in potentially
  multiple sub-populations, an hence in an object of class
  \code{\link{multipleFilterResult}}. Accordingly, \code{curv2Filters}
  can be used to split flow cytometry data sets.

}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{curv2Filters}.

}


\author{ Florian Hahne }


\seealso{

  \code{\link[flowCore:curv1Filter-class]{curv1Filter}},
  \code{\link{flowFrame}}, 
  \code{\link{flowSet}}, \code{\link[flowCore:filter-methods]{filter}}
  for evaluation 
  of \code{curv2Filters} and \code{\link{split}} for splitting of flow
  cytometry data sets based on that.
  
}

\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Create directly. Most likely from a command line
curv2Filter("FSC-H", "SSC-H", filterId="myCurv2Filter")

## To facilitate programmatic construction we also have the following
c2f <- curv2Filter(filterId="myCurv2Filter", x=list("FSC-H", "SSC-H"),
bwFac=2)
c2f <- curv2Filter(filterId="myCurv2Filter", x=c("FSC-H", "SSC-H"),
bwFac=2)

## Filtering using curv2Filter
fres <- filter(dat, c2f)
fres
summary(fres)
names(fres)

## The result of curv2 filtering are multiple sub-populations
## and we can split our data set accordingly
split(dat, fres)

## We can limit the splitting to one or several sub-populations
split(dat, fres, population="rest")
split(dat, fres, population=list(keep=c("area 2", "area 3")))


curv2Filter("FSC-H", "SSC-H", filterId="test filter")

}

\keyword{classes}
\keyword{methods}
