\name{curv1Filter-class}
\docType{class}
\alias{curv1Filter-class}
\alias{curv1Filter} 
\alias{show,curv1Filter-method}


\title{Class "curv1Filter"}


\description{
  
  Class and constructor for data-driven
  \code{\link[flowCore:filter-methods]{filter}} objects 
  that selects high-density regions in one dimension.

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
      describing the parameter used to filter the \code{flowFrame}. }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("curvFilter",
    ...)} or using the constructor \code{curv1Filter}.  Using the
    constructor is the recommended way of object instantiation:

}


\usage{

curv1Filter(x, bwFac=1.2, gridsize=rep(401, 2),
filterId="defaultCurv1Filter")

}


\arguments{

  \item{x}{ Character giving the name of the measurement parameter on
    which the filter is supposed to work on. This can also be a list
    containing a single character scalar for programmatic access. }
    
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    slot of this filter. The object can later be identified by this
    name. }
  
  \item{bwFac, gridsize}{ Numerics of length 1 and 2, respectively, used
    to set the \code{bwFac} and \code{gridsize} slots of the object. }
  
}


\value{

  Returns a \code{curv1Filter} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"curv1Filter")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "curv1Filter")}: Print
      information about the filter. }
    
  }
}


\details{

  Areas of high local density in one dimensions are identified by
  detecting significant curvature regions. See \cite{Duong, T. and
  Cowling, A. and Koch, I. and Wand, M.P., Computational Statistics and
  Data Analysis 52/9, 2008} for details. The constructor
  \code{curv1Filter} is a convenience function for object
  instantiation. Evaluating a \code{curv1Filter} results in potentially
  multiple sub-populations, an hence in an object of class
  \code{\link{multipleFilterResult}}. Accordingly, \code{curv1Filters}
  can be used to split flow cytometry data sets. 

}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{curv1Filters}.

}


\author{ Florian Hahne }


\seealso{

  \code{\link[flowCore:curv2Filter-class]{curv2Filter}},
  \code{\link{flowFrame}}, \code{\link{flowSet}},
  \code{\link[flowCore:filter-methods]{filter}} for evaluation 
  of \code{curv1Filters} and \code{\link{split}} for splitting of flow
  cytometry data sets based on that.
  
}


\examples{


## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Create directly. Most likely from a command line
curv1Filter("FSC-H", filterId="myCurv1Filter", bwFac=2)

## To facilitate programmatic construction we also have the following
c1f <- curv1Filter(filterId="myCurv1Filter", x=list("FSC-H"), bwFac=2)

## Filtering using curv1Filter
fres <- filter(dat, c1f)
fres
summary(fres)
names(fres)

## The result of curv1 filtering are multiple sub-populations
## and we can split our data set accordingly
split(dat, fres)

## We can limit the splitting to one or several sub-populations
split(dat, fres, population="rest")
split(dat, fres, population=list(keep=c("peak 2", "peak 3")))

}


\keyword{classes}
\keyword{methods}
