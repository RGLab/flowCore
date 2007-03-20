\name{ellipsoidGate-class}
\docType{class}
\alias{ellipsoidGate-class}
\alias{ellipsoidGate}
   
\title{Define ellipsoidal filter boundaries}
\description{
  A convenience method to facilitate the construction of a ellipsoid
  filter \code{\link{filter}} objects. Ellipsoid gates in n dimensions
  (n >= 2) are specified by their two foci (n coordinates each) and a distance
(not negative number).
}
\usage{
ellipsoidGate(filterId = "ellipsoidGate", .gate, distance, \dots)
}

\arguments{
  \item{filterId}{ An optional parameter that sets the \code{filterId} of this gate. }
  \item{.gate}{ A definition of the gate or foci. This can be either a list or matrix, described below. }
  \item{distance}{Positive numeric value}
  \item{\dots}{ You can also directly describe a gate without wrapping it in a list, as described below. }
  
}
\details{
	This function is designed to be useful in both direct and programmatic usage. In the first case, simply describe 
	the limits on the desired parameter directly as in the first example. To use this function programmatically, you
	may either construct a list using the same format that the direct description or you may construct a matrix with
	\code{n} columns and \code{2} rows. The first row corresponds to the minimal value for each parameter while the 
	second row corresponds to the maximal value for each parameter. The names of the parameters are taken from the 
	column names as in the third example.
}
\value{
	Returns a \code{\link{ellipsoidGate}} object for use in filtering \code{\link{flowFrame}}s or 
	other flow cytometry objects.
}
\author{ B. Ellis N. Le Meur}
\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }
\examples{
##Create directly. Most likely from a command line
##ellipsoidGate(filterId= "ellipsoidGateI", "FSC-H"= c(300, 200),
##"SSC-H"= c(500, 200), distance=400)

#To facilitate programmatic construction we also have the following
##ellipsoidGate(filterId = "ellipsoidGateII", list("FSC-H"= c(300, 200),
##"SSC-H" = c(300, 200)), distance=400)

##m <- matrix(c(300, 200, 500, 200),ncol = 2)
##colnames(m) <- c("FSC-H", "SSC-H")
##ellipsoidGate(filterId= "rectangleGate", m, distance=400)

}
\keyword{methods}

