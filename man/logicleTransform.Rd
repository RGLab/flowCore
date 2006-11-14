\name{logicleTransform}

\title{Create a 'Logicle' transformation object}
\description{Creates a \code{\link[flowCore]{biexponentialTransformation}} object with
parameters selected to closely match those described in the 'Logicle' transformation of
Parks et al. 
}

\usage{
logicleTransformation(w=0,r=262144,d=5,...)
}
\arguments{
  \item{w}{The width of the negative region (in decades). This defaults to 0.}
  \item{r}{The range of the ADC converter. This defaults to an 18-bit ADC with a range of 262144, appropriate for LSRII data.}
  \item{d}{The number of decades represented by the range of the ADC. This defaults to 5, the appropriate range for LSRII data.}
  \item{...}{Additional parameters for the \code{\link[flowCore]{biexponentialTransformation}} class.}
}

\details{
}

\value{Object of class \code{\link[flowCore]{transformation}.}
}  

\seealso{\code{\link[flowCore::transformation-class]{transformation-class}}
}

\references{
	Parks, David R., Roederer, Mario, Moore, Wayne A. (2006) \emph{Cytometry Part A} \bold{69A}, 541--551
}
\author{B. Ellis}
\examples{
  ##samp <- read.FCS(system.file("data","0877408774.B08", package="flowCore"))
  ##logicle <- logicleTransform(transformationId="Logicle-transformation",parameters="FSC-H",w=1)
  ##dataTransform <- applyTransformation(logicle,samp)
}
\keyword{methods}  

