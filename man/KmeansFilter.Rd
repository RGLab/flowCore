\name{kmeansFilter-class}
\docType{class}
\alias{kmeansFilter}
\alias{kmeansFilter-class}
\alias{length,kmeansFilter-method}
\alias{show,kmeansFilter-method}

\title{One dimensional k-means filters}

\description{

  A filter that performs one-dimensional k-means (Lloyd-Max) clustering
  on a parameter
}

\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("kmeansFilter",
    ...)} or using the constructor \code{kmeansFilter}.
  Using the constructor is the recommended way:

}

\usage{

  kmeansFilter(\dots, filterId = "kmeans")

}

\arguments{
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this gate. }
  
  \item{\dots}{ This gate takes a parameter name and a list of
    populations. They can be given as a character vector or as a named
    list.}
  
}

\section{Slots}{ 
  \describe{

    \item{\code{populations}:}{Object of class \code{"character"}. The
      names of the populations used for the \code{k} clusters to be
      identified by the filter. }
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameter used to filter the \code{flowFrame}. }
    
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
	"kmeansFilter")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user. }
    
    \item{show}{\code{signature(object = "kmeansFilter")}: Print
      information about the filter. }
    
  }
}

\details{

  The one-dimensional k-means filter is a multiple population filter
  capable of operating on a single parameter. It takes a parameter
  argument associated with two or more populations. Populations are
  considered to be ordered such that the population with the smallest
  mean intensity will be the first population in the list and the
  population with the highest mean intensity will be the last population
  listed.

}

\value{

  Returns a \code{\link{kmeansFilter}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}
      
\author{ B. Ellis, N. LeMeur, F.Hahne}

\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }

\examples{

a.1 <- kmeansFilter("FSC-H"=c("Pop1","Pop2","Pop3"), filterId="test filter")
summary(filter(read.FCS(system.file("extdata","0877408774.B08",
package="flowCore")),a.1))
}

\keyword{methods}
\keyword{classes}

