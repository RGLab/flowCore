\name{norm2Filter-class}
\docType{class}
\alias{norm2Filter-class}
\alias{norm2Filter}
\alias{show,norm2Filter-method}


\title{Class "norm2Filter"}


\description{

  Class and constructors for a \code{\link{filter}} that fits a
  bivariate normal distribution to a data set of paired values and
  selects data points according to their standard deviation from the
  fitted distribution.

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
    
    \item{\code{method}:}{One of \code{covMcd} or \code{cov.rob}
      defining method used for computation of covariance matrix.}
    
    \item{\code{scale.factor}:}{Numeric vector giving factor of standard
      deviations used for data selection (all points within
      \code{scalefac} standard deviations are selected).}

    \item{\code{n}:}{Object of class \code{"numeric"}, the number of
      events used to compute the covariance matrix of the bivariate
      distribution.}
  
    \item{\code{filterId}:}{Object of class \code{"character"}
      referencing the filter. }

    \item{\code{parameters}:}{Object of class \code{"ANY"} describing
      the parameters used to filter the \code{\link{flowFrame}} or
      \code{\link{flowSet}}.}

  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("norm2Filter",
    ...)} or using the constructor \code{norm2Filter}. The constructor
  is the recommended way of object instantiation:

}


\usage{

norm2Filter(x, y, method="covMcd", scale.factor=1, n=50000,
filterId="defaultNorm2Filter")

}


\arguments{

  \item{x,y}{ Characters giving the names of the measurement parameter
    on which the filter is supposed to work on. \code{y} can be missing
    in which case \code{x} is expected to be a character vector of
    length 2 or a list of characters. }
    
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    slot of this filter. The object can later be identified by this
    name.}
  
  \item{scale.factor, n}{ Numerics of length 1, used to set the
    \code{scale.factor} and \code{n} slots of the object. }

  \item{method}{ Character in \code{covMcd} or \code{cov.rob}, used to
    set the \code{method} slot of the object. }
  
}


\value{

  Returns a \code{\link{norm2Filter}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"norm2Filter")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "norm2Filter")}: Print
      information about the filter. }
    
  }
}


\details{

  The filter fits a bivariate normal distribution to the data and
  selects all events within the Mahalanobis distance multiplied by the
  \code{scale.factor} argument. The constructor \code{norm2Filter} is a
  conveniance function for object instantiation. Evaluating a
  \code{curv2Filter} results in an object of class
  \code{\link{logicalFilterResult}}. Accordingly, \code{norm2Filters}
  can be used to subset and to split flow cytometry data sets.

}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{norm2Filters}.

}


\author{F. Hahne}

\seealso{

  \code{\link[MASS::cov.rob]{cov.rob}},
  \code{\link[rrcov::covMcd]{covMcd}}, \code{\link[flowCore]{filter}}
  for evaluation of \code{norm2Filters} and \code{\link{split}} and
  \code{\link{Subset}}for splitting and subsetting of flow cytometry
  data sets based on that.

}

\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Create directly. Most likely from a command line
norm2Filter("FSC-H", "SSC-H", filterId="myCurv2Filter")

## To facilitate programmatic construction we also have the following
n2f <- norm2Filter(filterId="myNorm2Filter", x=list("FSC-H", "SSC-H"),
scale.factor=2)
n2f <- norm2Filter(filterId="myNorm2Filter", x=c("FSC-H", "SSC-H"),
scale.factor=2)

## Filtering using norm2Filter
fres <- filter(dat, n2f)
fres
summary(fres)

## The result of norm2 filtering is a logical subset
Subset(dat, fres)

## We can also split, in which case we get those events in and those
## not in the gate as separate populations
split(dat, fres)

}

\keyword{classes}
\keyword{methods}
