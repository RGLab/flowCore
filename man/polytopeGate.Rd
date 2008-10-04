\name{polytopeGate-class}
\docType{class}
\alias{polytopeGate-class}
\alias{polytopeGate}
\alias{show,polytopeGate-method}

   
\title{Define filter boundaries}

\description{
  
Convenience methods to facilitate the construction of
\code{\link{filter}} objects

}

\usage{
polytopeGate(filterId = "Polytope Gate", .gate, \dots)

}

\arguments{
  
  \item{filterId}{ An optional parameter that sets the \code{filterId}
    of this gate. }
  
  \item{.gate}{ A definition of the gate. This can be either a list,
    vector or matrix, described below. }
  
  
  \item{\dots}{ You can also directly describe a gate without wrapping
    it in a list or matrix, as described below. }
  
} \details{

  These functions are designed to be useful in both direct and
  programmatic usage.
  
  For rectangle gate in n dimensions, if n=1 the gate correspond to a
  range gate. If n=2, the gate is a rectangle gate. To use this function
  programmatically, you may either construct a list or you may construct
  a matrix with \code{n} columns and \code{2} rows.  The first row
  corresponds to the minimal value for each parameter while the second
  row corresponds to the maximal value for each parameter.  The names of
  the parameters are taken from the column names as in the third
  example.

  Rectangle gate objects can also be multiplied together using the
  \code{*} operator, provided that both gate have orthogonal axes.

  For polygone gate, the boundaries are specified as vertices in 2
  dimensions, for polytope gate objects as vertices in n dimensions.

  For quadrant gates, the boundaries are specified as a named list or
  vector of length two.
  
}
\value{
  
  Returns a \code{\link{rectangleGate}} or \code{\link{polygonGate}}
  object for use in filtering \code{\link{flowFrame}}s or other flow
  cytometry objects.
  
}

\author{F.Hahne, B. Ellis N. Le Meur}

\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }

\examples{


}
\keyword{methods}

