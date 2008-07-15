\name{timeFilter-class}
\docType{class}
\alias{timeFilter-class}
\alias{timeFilter} 
\alias{timeFilter-class}
\alias{show,timeFilter-method}

\title{Class "timeFilter" }

\description{

  Define a filter that removes stretches of unusual data distribution
  within a single parameter over time. This can be used to correct for
  problems during data acquisition like air bubbles or clods.

}

\section{Objects from the Class}{

  Objects can be created by calls of the form
  \code{new("timeFilter", ...)} or using the constructor
  \code{timeFilter}. Using the constructor is the recommended way:

}


\usage{

timeFilter(\dots, filterId="time", bandwidth=0.75, binSize=NULL,
           timeParameter=NULL)

}

\arguments{
  
  \item{\dots}{ The names of the parameters on which the filter is
    supposed to work on. Names can either be given as individual
    arguments, or as a list or a character vector. }
    
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    slot of this gate. }
  
  \item{bandwidth, binSize}{ Numerics used to set the \code{bandwidth}
    and \code{binSize} slots of the object. }
  
  \item{timeParameter}{ Character used to set the \code{timeParameter}
    slot of the object.}
  
}

\section{Slots}{ 
  \describe{

    \item{\code{bandwidth}:}{Object of class \code{"numeric"}. The
      sensitivity of the filter, i.e., the amount of local variance of
      the signal we want to allow. }
    
    \item{\code{binSize}:}{Object of class \code{"numeric"}. The size
      of the bins used for the local variance and location
      estimation. If \code{NULL}, a reasonable default is used when
      evvaluating the filter. }

    \item{\code{timeParameter}:}{Object of class \code{"character"},
      used to define the time domain parameter. If \code{NULL}, the
      filter tries to guess the time domain from the  \code{flowFrame}.}
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameters used to filter the \code{flowFrame}. }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}

\section{Extends}{
  
  Class \code{"\linkS4class{parameterFilter}"}, directly.
  Class \code{"\linkS4class{concreteFilter}"}, by class
  \code{parameterFilter}, distance 2.
  Class \code{"\linkS4class{filter}"}, by class \code{parameterFilter},
  distance3.
  
}

\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"timeFilter")}: The workhorse used to evaluate the filter on
      data. This is usually not called directly by the user. }
      
    \item{show}{\code{signature(object = "timeFilter")}: Print
      information about the filter. }
    
  }
}


\author{ Florian Hahne }


\seealso{
  \code{\link[flowCore]{filter}}
}

\examples{
timeFilter("SSC-H", bandwidth=1)

\keyword{classes}
