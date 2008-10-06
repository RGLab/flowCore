\name{kmeansFilter-class}
\docType{class}
\alias{kmeansFilter}
\alias{kmeansFilter-class}
\alias{length,kmeansFilter-method}
\alias{show,kmeansFilter-method}


\title{Class "kmeansFilter"}


\description{

  A filter that performs one-dimensional k-means (Lloyd-Max) clustering
  on a single flow parameter. }


\section{Extends}{
  
  Class \code{\linkS4class{parameterFilter}}, directly.

  Class \code{\linkS4class{concreteFilter}}, by class
  \code{parameterFilter}, distance 2.

  Class \code{\linkS4class{filter}}, by class \code{parameterFilter},
  distance3.

}


\section{Slots}{ 
  \describe{

    \item{\code{populations}:}{Object of class \code{character}. The
      names of the \code{k} populations (or clusters) that will be
      created by the \code{kmeansFilter}. These names will later be used
      for the respective subpopulations in \code{\link{split}}
      operations and for the summary of the \code{\link{filterResult}}. }
    
    \item{\code{parameters}:}{Object of class \code{\link{parameters}},
      defining a single parameter for which the data in the
      \code{\linkS4class{flowFrame}} is to be clustered. This may also
      be a \code{\link[transform-class]{transformation}} object. }
    
    \item{\code{filterId}:}{Object of class \code{character}, an
      identifier or name to reference the \code{kmeansFilter} object
      later on. }
    
  }
}


\section{Objects from the Class}{
  
  Like all other \code{\linkS4class{filter}} objects in \code{flowCore},
  \code{kmeansFilter} objects should be instantiated through their
  constructor \code{kmeansFilter()}. See the \code{Usage} section for
  details.
  
}


\usage{

kmeansFilter(\dots, filterId="defaultKmeansFilter")

}

\arguments{
  
  \item{\dots}{ \code{kmeansFilter} are defined by a single flow
    parameter and an associated list of \code{k} population names. They
    can be given as a character vector via a named argument, or as a
    list with a single named argument. In both cases the name will be
    used as the flow parameter and the content of the list or of the
    argument will be used as population names, after coercing to
    character. For example

    \code{kmeansFilter(FSC=c("a", "b", "c"))}

    or

    \code{kmeansFilter(list(SSC=1:3))}

    If the parameter is not fully realized, but instead is the result of a
    \code{\link[transform-class]{transformation}} operation, two
    arguments need to be passed to the constructor: the first one being
    the \code{\link[transform-class]{transform}} object and the second
    being a vector of population names which can be coerced to a
    character. For example

    \code{kmeansFilter(tf, c("D", "E"))}
    
    
  }
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of the object. The filter can later be identified by this name. }
 
}


\value{

  Returns a \code{kmeansFilter} object for use in filtering
  \code{\link[flowFrame-class]{flowFrames}} or other flow cytometry
  objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"kmeansFilter")}: The workhorse used to evaluate the filter on
      data.

      \emph{Usage:}
      
      This is usually not called directly by the user, but internally by
      the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "kmeansFilter")}: Print
      information about the filter.

      \emph{Usage:}

      The method is called automatically whenever the object is printed
      on the screen. }
    
  }
}


\details{

  The one-dimensional k-means filter is a multiple population filter
  capable of operating on a single flow parameter. It takes a parameter
  argument associated with two or more populations and results in the
  generation of an object of class \code{\link{multipleFilterResult}}.
  Populations are considered to be ordered such that the population with
  the smallest mean intensity will be the first population in the list
  and the population with the highest mean intensity will be the last
  population listed.

}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{kmeansFilters}.

}


\author{ F. Hahne, B. Ellis, N. LeMeur}


\seealso{

  \code{\link{flowFrame}}, \code{\link{flowSet}}, \code{\link{filter}}
  for evaluation of \code{kmeansFilters} and \code{\link{split}} for
  splitting of flow cytometry data sets based on the result of the
  filtering operation.

}


\examples{
## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Create the filter
kf <- kmeansFilter("FSC-H"=c("Pop1","Pop2","Pop3"), filterId="myKmFilter")

## Filtering using kmeansFilters
fres <- filter(dat, kf)
fres
summary(fres)
names(fres)

## The result of quadGate filtering are multiple sub-populations
## and we can split our data set accordingly
split(dat, fres)

## We can limit the splitting to one or several sub-populations
split(dat, fres, population="Pop1")
split(dat, fres, population=list(keep=c("Pop1","Pop2")))

}

\keyword{methods}
\keyword{classes}

