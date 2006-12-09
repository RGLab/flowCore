\name{rectangleGate}
\alias{rectangleGate}


\title{Define filter boundaries}
\description{
  A convenience method to facilitate the construction of \code{\link{filter}} objects
}
\usage{
rectangleGate(filterId = "rectangleGate", .gate, \dots)
polygonGate(filterId = "polygonGate", .gate, \dots)
polytope(filterId = "polygonGate", .gate, \dots)
ellipsoidGate(filterId = "polygonGate", .gate, distance, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{ An optional parameter that sets the \code{filterId} of this
    gate. }
  \item{.gate}{ A definition of the gate. This can be either a list or matrix, described below. }
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
	Returns a \code{\link{rectangleGate}} object for use in filtering \code{\link{flowFrame}}s or 
	other flow cytometry objects.
}
\author{ B. Ellis }
\seealso{ \code{\link{flowFrame}}, \code{\link{filter}} }
\examples{
#Create directly. Most likely from a command line
rectangleGate(filterId="foo","FSC-H"=c(1,50),"SSC-H"=c(50,800))

#To facilitate programmatic construction we also have the following
rectangleGate(filterId="foo",list("FSC-H"=c(1,50),"SSC-H"=c(50,800)))

m = matrix(c(1,50,50,800),ncol=2)
colnames(m) = c("FSC-H","SSC-H")
rectangleGate(filterId="foo",m)

}
\keyword{ misc }

