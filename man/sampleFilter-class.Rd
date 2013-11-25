\name{sampleFilter-class}
\docType{class}
\alias{sampleFilter-class}
\alias{sampleFilter}
\alias{show,sampleFilter-method}


\title{Class "sampleFilter"}


\description{
  
  This non-parameter filter selects a number of events from the primary
  \code{\link{flowFrame}}.

}


\section{Extends}{
  
  Class \code{"\linkS4class{concreteFilter}"}, directly.

  Class \code{"\linkS4class{filter}"}, by class
  \code{concreteFilter}, distance 2.

}


\section{Slots}{
  \describe{
    
    \item{\code{size}:}{Object of class \code{"numeric"}. Then number of
      events that are to be selected.}
    
    \item{\code{filterId}:}{A character vector that identifies this
      \code{filter}.}
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("sampleFilter",
    ...)} or using the constructor \code{sampleFilter}. The latter is the
    recommended way of object instantiation:

}


\usage{

sampleFilter(size, filterId="defaultSampleFilter")

}

\arguments{

  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this \code{\link{filter}}. The object can later be identified by
    this name.}
  
  \item{size}{The number of events to select.}
  
}


\value{
  
  Returns a \code{sampleFilter} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"sampleFilter")}: The workhorse used to evaluate the gate on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "sampleFilter")}: Print
      information about the gate. }
    
  }
}


\details{
  
  Selects a number of events without replacement from a
  \code{flowFrame}.
  
}


\author{B. Ellis, F.Hahne}


\seealso{

  \code{\link{flowFrame}}, \code{\link{filter}} for evaluation of
  \code{sampleFilters} and \code{\link{split}} and
  \code{\link{Subset}}for splitting and subsetting of flow cytometry
  data sets based on that.
  
}

\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

#Create the filter
sf <- sampleFilter(filterId="mySampleFilter", size=500)
sf

## Filtering using sampeFilters
fres <- filter(dat, sf)
fres
summary(fres)

## The result of sample filtering is a logical subset
Subset(dat, fres)

## We can also split, in which case we get those events in and those
## not in the gate as separate populations
split(dat, fres)

}


\keyword{methods}
\keyword{classes}
