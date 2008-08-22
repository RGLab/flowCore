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
polytopeGate(filterId = "polytopeGate", .gate, \dots)

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
data(GvHD)

#Create directly. Most likely from a command line
rectangleGate(filterId="rectangleI", "FSC-H"=c(1, 50), "SSC-H"=c(50, 800))
quadGate(filterId="quadI", "FSC-H"=100, "SSC-H"=400)

#To facilitate programmatic construction we also have the following
rectangleGate(filterId="rectangleII", list("FSC-H"=c(1, 50),
"SSC-H"=c(50, 800)))
quadGate(filterId="quadII", list("FSC-H"=100, "SSC-H"=400))

m <- matrix(c(50,100,50,700),ncol=2)
colnames(m) <- c("FSC-H","SSC-H")
rect <- rectangleGate(filterId="rectangleGate", m)
pop1 <- filter(GvHD[1],rect)
sub1 <- Subset(GvHD[1], pop1)

## Multiply rectangle gates
rect2 <- rectangleGate(filterId="FSC-", "FSC-H"=c(-Inf, 50)) *
rectangleGate(filterId="SSC+", "SSC-H"=c(50, Inf))
pop2 <- filter(GvHD[1],rect2)
sub2 <- Subset(GvHD[1], pop2)


## PolygonGate
sqrcut <- matrix(c(300,300,600,600,50,300,300,50),ncol=2,nrow=4)
colnames(sqrcut) <- c("FSC-H","SSC-H")
pg <- polygonGate(filterId="nonDebris", boundaries= sqrcut)
pg
pop3 <- filter(GvHD[1], pg)
sub3 <- Subset(GvHD[1], pop3)

## QuadGate
qg <- quadGate(filterId="quad", "FSC-H"=100, "SSC-H"=400)
pop4 <- filter(GvHD[1], qg)
sub4 <- split(GvHD[1], pop4)


}
\keyword{methods}

