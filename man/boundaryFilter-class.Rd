\name{boundaryFilter-class}
\docType{class}
\alias{boundaryFilter-class}
\alias{boundaryFilter} 
\alias{show,boundaryFilter-method}


\title{Class "boundaryFilter"}


\description{
  
  Class and constructor for data-driven \code{\link{filter}} objects
  that discard margin events.

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

    \item{\code{tolerance}:}{Object of class \code{"numeric"}. The
      machine tolerance used to decide whether an event is on the
      measurement boundary. Essentially, this is done by evaluating
      \code{x>minRange+tolerance & x<maxRange-tolerance}. }

     \item{\code{side}:}{Object of class \code{"character"}. The margin
       on which to evaluate the filter. Either \code{upper} for the
       upper margin or \code{lower} for the lower margin or \code{both}
       for both margins. }
      
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("boundaryFilter",
    ...)} or using the constructor \code{boundaryFilter}.  Using the
    constructor is the recommended way of object instantiation:

}


\usage{

boundaryFilter(x, tolerance=.Machine$double.eps, side=c("both", "lower",
"upper"), filterId="defaultBoundaryFilter")

}


\arguments{

  \item{x}{ Character giving the name(s) of the measurement parameter(s) on
    which the filter is supposed to work. Note that all events on the
    margins of ay of the channels provided by \code{x} will be
    discarded, which is often not desired. Such events may not
    convey much information in the particular channel on which their
    value falls on the margin, however they may well be informative in
    other channels. }

  \item{tolerance}{ Numeric vector, used to set the \code{tolerance}
    slot of the object. Can be set separately for each element in
    \code{x}. R's recycling rules apply. }

   \item{side}{ Character vector, used to set the \code{side}
    slot of the object.  Can be set separately for each element in
    \code{x}. R's recycling rules apply.}
    
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    slot of this filter. The object can later be identified by this
    name. }
  
 
  
}


\value{

  Returns a \code{boundaryFilter} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"boundaryFilter")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "boundaryFilter")}: Print
      information about the filter. }
    
  }
}


\details{

  Flow cytomtery instruments usually operate on a given data range, and
  the limits of this range are stored as keywords in the FSC
  files. Depending on the amplification settings and the dynamic range
  of the measured signal, values can occur that are outside of the
  measurement range, and most instruments will simply pile those values
  at the minimum or maximum range limit. The \code{boundaryFilter}
  removes these values, either for a single parameter, or for a
  combination of parameters. Note that it is often desirable to treat
  boundary events on a per-parameter basis, since their values might be
  uninformative for one particular channel, but still be useful in all of
  the other channels.

  The constructor \code{boundaryFilter} is a convenience function for
  object instantiation. Evaluating a \code{boundaryFilter} results in a
  single sub-populations, an hence in an object of class
  \code{\link{filterResult}}.

}



\author{ Florian Hahne }


\seealso{

 \code{\link{flowFrame}}, \code{\link{flowSet}},
  \code{\link[flowCore:filter-methods]{filter}} for evaluation of
  \code{boundaryFilters} and \code{\link{Subset}} for subsetting of flow
  cytometry data sets based on that.
  
}


\examples{


## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Create directly. Most likely from a command line
boundaryFilter("FSC-H", filterId="myBoundaryFilter")

## To facilitate programmatic construction we also have the following
bf <- boundaryFilter(filterId="myBoundaryFilter", x=c("FSC-H"))

## Filtering using boundaryFilter
fres <- filter(dat, bf)
fres
summary(fres)

## We can subset the data with the result from the filtering operation.
Subset(dat, fres)

## A boundaryFilter on the lower margins of several channels
bf2 <- boundaryFilter(x=c("FSC-H", "SSC-H"), side="lower")

}


\keyword{classes}
\keyword{methods}
