\name{quadGate-class}
\docType{class}
\alias{quadGate-class}
\alias{quadGate}
\alias{show,quadGate-method}


\title{Class "quadGate"}


\description{

  Class and constructors for quadrant-type \code{\link{filter}} objects.

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

    \item{\code{boundary}:}{Object of class \code{"numeric"}, length
      2. The boundaries of the quadrant regions.}
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      describing the parameter used to filter the \code{flowFrame}. }
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the gate.}
    
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("quadGate",
    ...)} or using the constructor \code{quadGate}. The latter is the
    recommended way of object instantiation:

}


\usage{

quadGate(filterId = "Quadrant Gate", .gate, \dots)

}

\arguments{
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this \code{\link{filter}}. The object can later be identified by
    this name.}
  
  \item{.gate}{ A definition of the gate for programmatic access. This
    can be either a named list or a named numeric vector, as described
    below. }
  
  \item{\dots}{ The parameters of \code{quadGates} can also be directly
    described using named function arguments, as described below. }
  
}


\value{
  
  Returns a \code{quadGate} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"quadGate")}: The workhorse used to evaluate the gate on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "quadGate")}: Print
      information about the gate. }
    
  }
}


\details{

  \code{quadGates} are defined by two parameters, which specify a
  separation of a two-dimensional parameter space into four
  quadrants. The \code{quadGate} function is designed to be useful in
  both direct and programmatic usage:

  For the interactive use, these parameters can be given as additional
  named function arguments, where the names correspond to valid
  parameter names in a \code{\link{flowFrame}} or
  \code{\link{flowSet}}. For a more programmatic approach, a named list
  or numeric vector of the gate boundaries can be passed on to the
  function as argument \code{.gate}.

  Evaluating a \code{quadGate} results in four sub-populations, and hence
  in an object of class \code{\link{multipleFilterResult}}. Accordingly,
  \code{quadGates} can be used to split flow cytometry data sets.
  
}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for plotting of
  \code{quadGates}.

}


\author{F.Hahne, B. Ellis N. Le Meur}


\seealso{

  \code{\link{flowFrame}}, \code{\link{flowSet}}, \code{\link{filter}}
  for evaluation of \code{quadGates} and \code{\link{split}} for
  splitting of flow cytometry data sets based on that.

}


\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

## Create directly. Most likely from a command line
quadGate(filterId="myQuadGate1", "FSC-H"=100, "SSC-H"=400)

## To facilitate programmatic construction we also have the following
quadGate(filterId="myQuadGate2", list("FSC-H"=100, "SSC-H"=400))
## FIXME: Do we want this?
##quadGate(filterId="myQuadGate3", .gate=c("FSC-H"=100, "SSC-H"=400))

## Filtering using quadGates
qg <- quadGate(filterId="quad", "FSC-H"=600, "SSC-H"=400)
fres <- filter(dat, qg)
fres
summary(fres)
names(fres)

## The result of quadGate filtering are multiple sub-populations
## and we can split our data set accordingly
split(dat, fres)

## We can limit the splitting to one or several sub-populations
split(dat, fres, population="FSC-H-SSC-H-")
split(dat, fres, population=list(keep=c("FSC-H-SSC-H-",
"FSC-H-SSC-H+")))

}


\keyword{classes}
\keyword{methods}

