\name{ellipsoidGate-class}
\docType{class}
\alias{ellipsoidGate-class}
\alias{ellipsoidGate}
\alias{show,ellipsoidGate-method}

\title{Class "ellipsoidGate"}

\description{

  Class and constructor for n-dimensional ellipsoidal
  \code{\link{filter}} objects.
  
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

    \item{\code{mean}:}{Objects of class \code{"numeric"}. Vector giving
      the location of the center of the ellipse in n dimensions. }

    \item{\code{cov}:}{Objects of class \code{"matrix"}. The covariance
      matrix defining the shape of the ellipse. }

    \item{\code{distance}:}{Objects of class \code{"numeric"}. The
      Mahalanobis distance defining the size of the ellipse. }
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameter used to filter the \code{flowFrame}. }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("ellipsoidGate",
    ...)} or by using the constructor \code{ellipsoidGate}.  Using the
    constructor is the recommended way of object instantiation:

}


\usage{

ellipsoidGate(\dots, .gate, mean, distance=1, filterId="defaultEllipsoidGate")

}


\arguments{
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this gate.}
  
  \item{.gate}{ A definition of the gate via a covariance matrix. }
  
  \item{mean}{ Numeric vector of equal length as dimensions in
    \code{.gate}.}

  \item{distance}{ Numeric scalar giving the Mahalanobis distance
      defining the size of the ellipse. This mostly exists for
      compliance reasons to the gatingML standard as \code{mean} and
      \code{gate} should already uniquely define the
      ellipse. Essentially, \code{distance} is merely a factor that gets
      applied to the values in the covariance matrix.}
  
  \item{\dots}{ You can also directly describe the covariance matrix
    through named arguments, as described below. }
  
}


\value{
	
  Returns a \code{\link{ellipsoidGate}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
	
}


\details{

  A convenience method to facilitate the construction of a ellipsoid
  \code{\link{filter}} objects. Ellipsoid gates in n dimensions (n >= 2)
  are specified by a a covarinace matrix and a vector of mean values
  giving the center of the ellipse.
  
  This function is designed to be useful in both direct and programmatic
  usage. In the first case, simply describe the covariance matrix
  through named arguments. To use this function programmatically, you
  may pass a covarince matrix and a mean vector directly, in which case
  the parameter names are the colnames of the matrix.

}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"ellipsoidGate")}: The workhorse used to evaluate the filter on
      data. This is usually not called directly by the user, but
      internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "ellipsoidGate")}: Print
      information about the filter. }
    
  }
}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{ellipsoidGates}.

}


\author{F.Hahne, B. Ellis, N. LeMeur}

\seealso{

  \code{\link{flowFrame}}, \code{\link{polygonGate}},
  \code{\link{rectangleGate}}, \code{\link{polytopeGate}},
  \code{\link{filter}} for evaluation of \code{rectangleGates} and
  \code{\link{split}} and \code{\link{Subset}}for splitting and
  subsetting of flow cytometry data sets based on that.
  
}

\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Defining the gate
cov <- matrix(c(6879, 3612, 3612, 5215), ncol=2,
dimnames=list(c("FSC-H", "SSC-H"), c("FSC-H", "SSC-H")))
mean <- c("FSC-H"=430, "SSC-H"=175)
eg <- ellipsoidGate(filterId= "myEllipsoidGate", .gate=cov, mean=mean)

## Filtering using ellipsoidGates
fres <- filter(dat, eg)
fres
summary(fres)

## The result of ellipsoid filtering is a logical subset
Subset(dat, fres)

## We can also split, in which case we get those events in and those
## not in the gate as separate populations
split(dat, fres)

##ellipsoidGate can be converted to polygonGate by interpolation
pg <- as(eg, "polygonGate")
pg


}
\keyword{methods}

