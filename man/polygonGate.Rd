\name{polygonGate-class}
\docType{class}
\alias{polygonGate-class}
\alias{polygonGate}
\alias{show,polygonGate-method}


\title{Class "polygonGate"}


\description{

  Class and constructor for 2-dimensional polygonal \code{\link{filter}}
  objects.
  
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

    \item{\code{boundaries}:}{Object of class \code{"matrix"}. The
      vertices of the polygon in two dimensions. There need to be at
      least 3 vertices specified for a valid polygon.}
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameter used to filter the \code{flowFrame}. }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("polygonGate",
    ...)} or by using the constructor \code{polygonGate}. Using the
    constructor is the recommended way of object instantiation:

}


\usage{
polygonGate(filterId = "polygonGate", boundaries, \dots)

}


\arguments{
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this gate. }
  
  \item{boundaries}{ A definition of the gate. This can be either a list
    or a named matrix as described below. }
  
  
  \item{\dots}{ You can also directly describe a gate without wrapping
    it in a list or matrix, as described below. }
  
}


\value{
  
  Returns a \code{\link{polygonGate}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"polygonGate")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "polygonGate")}: Print
      information about the filter. }
    
  }
}


\details{

  Polygons are specified by the coordinates of their vertices in two
  dimensions. The constructor is designed to be useful in both direct
  and programmatic usage. It takes either a list or a named matrix with
  \code{2} columns and at least \code{3} rows containing these
  coordinates. Alternatively, vertices can be given as named arguments,
  in which case the function tries to convert the values into a matrix.
  
}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{polygonGates}.

}


\author{F.Hahne, B. Ellis N. Le Meur}


\seealso{

  \code{\link{flowFrame}}, \code{\link{rectangleGate}},
  \code{\link{ellipsoidGate}}, \code{\link{polytopeGate}},
  \code{\link{filter}} for evaluation of \code{rectangleGates} and
  \code{\link{split}} and \code{\link{Subset}}for splitting and
  subsetting of flow cytometry data sets based on that.
  
}

\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Defining the gate
sqrcut <- matrix(c(300,300,600,600,50,300,300,50),ncol=2,nrow=4)
colnames(sqrcut) <- c("FSC-H","SSC-H")
pg <- polygonGate(filterId="nonDebris", boundaries= sqrcut)
pg

## Filtering using polygonGates
fres <- filter(dat, pg)
fres
summary(fres)

## The result of polygon filtering is a logical subset
Subset(dat, fres)

## We can also split, in which case we get those events in and those
## not in the gate as separate populations
split(dat, fres)

}
\keyword{methods}

