\name{rectangleGate-class}
\docType{class}
\alias{rectangleGate-class}
\alias{rectangleGate}
\alias{polygonGate-class}
\alias{polygonGate}
\alias{polytopeGate}
\alias{polytopeGate-class}
\alias{polytopeGate-class}
\alias{summary,rectangleGate-method}
\alias{show,rectangleGate-method}
\alias{show,polygonGate-method}
   
\title{Define filter boundaries}
\description{
  A convenience method to facilitate the construction of \code{\link{filter}} objects
}
\usage{
rectangleGate( filterId = "rectangleGate", .gate, \dots)
polygonGate(filterId = "polygonGate", boundaries, \dots)
polytopeGate(filterId = "polytopeGate", .gate, \dots)
}

\arguments{
  \item{filterId}{ An optional parameter that sets the \code{filterId} of this gate. }
  \item{.gate}{ A definition of the gate. This can be either a list or matrix, described below. }
  \item{boundaries}{ A definition of the gate.}
  \item{\dots}{ You can also directly describe a gate without wrapping it in a list, as described below. }
  
}
\details{
  This function is designed to be useful in both direct and programmatic usage.
  
  For rectangle gate in n dimensions, if n=1 the gate correspond to a
  range gate. If n=2, the gate is a rectangle gate. To use this function
  programmatically, you may either construct a list
  or you may construct a matrix with  \code{n} columns and \code{2}
  rows.   The first row corresponds to the minimal value for each parameter
  while the  second row corresponds to the maximal value for each
  parameter.  The names of the parameters are taken from the 
  column names as in the third example.

  For polygone gate, the boundaries are specified as vertices in 2 dimensions.
  
 
}
\value{
	Returns a \code{\link{rectangleGate}} or
	\code{\link{polygonGate}} object for use in filtering \code{\link{flowFrame}}s or 
	other flow cytometry objects.
}
\author{ B. Ellis N. Le Meur}
\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }
\examples{
data(GvHD)

#Create directly. Most likely from a command line
rectangleGate(filterId="rectangleGateI","FSC-H"=c(1,50),"SSC-H"=c(50,800))

#To facilitate programmatic construction we also have the following
rectangleGate(filterId="rectangleGateII",list("FSC-H"=c(1,50),"SSC-H"=c(50,800)))

m <- matrix(c(50,100,50,700),ncol=2)
colnames(m) <- c("FSC-H","SSC-H")
rect <- rectangleGate(filterId="rectangleGate", m)
pop1 <- filter(GvHD,rect)
sub3 <- Subset(GvHD, pop1)


## PolygonGate
sqrcut <- matrix(c(300,300,600,600,50,300,300,50),ncol=2,nrow=4)
colnames(sqrcut) <- c("FSC-H","SSC-H")
polygonGate(filterId="nonDebris", boundaries= sqrcut)

}
\keyword{methods}

