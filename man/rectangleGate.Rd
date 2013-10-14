\name{rectangleGate-class}
\docType{class}
\alias{rectangleGate-class}
\alias{rectangleGate}
\alias{summary,rectangleGate-method}
\alias{show,rectangleGate-method}
\alias{[,rectangleGate,character-method}
\alias{[,rectangleGate,ANY-method}
\alias{*,rectangleGate,rectangleGate-method}
   
\title{Class "rectangleGate"}

\description{

  Class and constructor for n-dimensional rectangular
  \code{\linkS4class{filter}} objects.

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

    \item{\code{min,max}:}{Objects of class \code{"numeric"}. The
      minimum and maximum values of the n-dimensional rectangular
      region. }
    
    \item{\code{parameters}:}{Object of class \code{"character"},
      indicating the parameters for which the \code{rectangleGate} is
      defined.}
    
    \item{\code{filterId}:}{Object of class \code{"character"},
      referencing the filter.}
    
  }
}


\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("rectangleGate",
    ...)}, by using the constructor \code{rectangleGate} or by combining
    existing \code{rectangleGates} using the \code{*} method.  Using the
    constructor is the recommended way of object instantiation:

}


\usage{

rectangleGate(\dots, .gate, filterId="defaultRectangleGate")

}

\arguments{
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this gate. The object can later be identified by this name.}
  
  \item{.gate}{ A definition of the gate. This can be either a list, or
    a matrix, as described below. }
  
  \item{\dots}{ You can also directly provide the boundaries of a
    \code{rectangleGate} as additional named arguments, as described
    below. }
  
}


\value{
  
  Returns a \code{\link{rectangleGate}} object for use in filtering
  \code{\link{flowFrame}}s or other flow cytometry objects.
  
}


\section{Methods}{
  \describe{
    
    \item{\%in\%}{\code{signature(x = "flowFrame", table =
	"rectangleGate")}: The workhorse used to evaluate the filter on
	data. This is usually not called directly by the user, but
	internally by calls to the \code{\link{filter}} methods. }
    
    \item{show}{\code{signature(object = "rectangleGate")}: Print
      information about the filter. }

    \item{*}{\code{signature(e1 = "rectangleGate", e2 =
	"rectangleGate")}: combining two \code{rectangleGates} into one
      higher dimensional representation. }

    \item{\[}{\code{signature(x = "rectangleGate", i = "character")}:
	Subsetting of a \code{rectangleGate} by parameter name(s). This
	is essentially the inverse to \code{*}. }
    
  }
}


\details{
  
  This class describes a rectangular region in n dimensions, which is a
  Cartesian product of \code{n} orthogonal intervals in these
  dimensions. \code{n=1} corresponds to a range gate, \code{n=2} to a
  rectangle gate, \code{n=3} corresponds to a box region and \code{n>3}
  to a hyper-rectangular regions. Intervals may be open on one side, in
  which case the value for the boundary is supposed to be \code{Inf} or
  \code{-Inf}, respectively. \code{rectangleGates} are inclusive, that
  means that events on the boundaries are considered to be in the gate.
  
  The constructor is designed to be useful in both direct and
  programmatic usage. To use it programmatically, you may either
  construct a named list or you may construct a matrix with \code{n}
  columns and \code{2} rows.  The first row corresponds to the minimal
  value for each parameter while the second row corresponds to the
  maximal value for each parameter.  The names of the parameters are
  taken from the column names or from the list names,
  respectively. Alternatively, the boundaries of the
  \code{rectangleGate} can be given as additional named arguments, where
  each of these arguments should be a numeric vector of length \code{2};
  the function tries to collapse these boundary values into a matrix.

  Note that boundaries of \code{rectangleGates} where \code{min > max}
  are syntactically valid, however when evaluated they will always be
  empty.

  \code{rectangleGate} objects can also be multiplied using the \code{*}
  operator, provided that both gates have orthogonal axes. This results
  in higher-dimensional \code{rectangleGates}. The inverse operation of
  subsetting by parameter name(s) is also available.

  Evaluating a \code{rectangleGate} generates an object of class
  \code{\linkS4class{logicalFilterResult}}. Accordingly,
  \code{rectangleGates} can be used to subset and to split flow
  cytometry data sets.
  
}


\note{

  See the documentation in the
  \code{\link[flowViz:flowViz-package]{flowViz}} package for details on
  plotting of \code{rectangleGates}.

}


\author{F.Hahne, B. Ellis N. Le Meur}


\seealso{

  \code{\link{flowFrame}}, \code{\link{polygonGate}},
  \code{\link{ellipsoidGate}}, \code{\link{polytopeGate}},
  \code{\link{filter}} for evaluation of \code{rectangleGates} and
  \code{\link{split}} and \code{\link{Subset}}for splitting and
  subsetting of flow cytometry data sets based on that.
  
}


\examples{

## Loading example data
dat <- read.FCS(system.file("extdata","0877408774.B08",
package="flowCore"))

#Create directly. Most likely from a command line
rectangleGate(filterId="myRectGate", "FSC-H"=c(200, 600), "SSC-H"=c(0, 400))

#To facilitate programmatic construction we also have the following
rg <- rectangleGate(filterId="myRectGate", list("FSC-H"=c(200, 600),
"SSC-H"=c(0, 400)))
mat <- matrix(c(200, 600, 0, 400), ncol=2, dimnames=list(c("min", "max"),
c("FSC-H", "SSC-H")))
rg <- rectangleGate(filterId="myRectGate", .gate=mat)

## Filtering using rectangleGates
fres <- filter(dat, rg)
fres
summary(fres)

## The result of rectangle filtering is a logical subset
Subset(dat, fres)

## We can also split, in which case we get those events in and those
## not in the gate as separate populations
split(dat, fres)

## Multiply rectangle gates
rg1 <- rectangleGate(filterId="FSC-", "FSC-H"=c(-Inf, 50))
rg2 <- rectangleGate(filterId="SSC+", "SSC-H"=c(50, Inf))
rg1 * rg2

## Subset rectangle gates
rg["FSC-H"]

}


\keyword{methods}
\keyword{classes}
