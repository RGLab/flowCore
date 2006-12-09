\name{filter}
\alias{filter}
\alias{filter,flowFrame,filter}
\alias{filter,flowFrame,rectangleGate}
\alias{filter,flowFrame,polygonGate}
\alias{filter,flowFrame,norm2Filter}

\title{Filter on a data set}
\description{filter is a function for subsetting data}

\usage{
filter(flowObject,filter)
}
\arguments{
  \item{flowObject}{Object of class \code{\link[flowCore]{flowFrame}}.}
  \item{filter}{Object of class gate and slots parameters and boundaries}
  \item{parent}{Logical. Assign center of distribution.}
}

\details{
}

\value{Object of class \code{\link[flowCore]{filterResult}}
}  

\seealso{\code{\link[flowCore::gate-class]{gate-class}},
  \code{\link[flowCore::filterResult-class]{filterResult-class}}
}

\author{N. Le Meur}
\examples{
 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 rectGate <- rectangleGate(filterId="nonDebris","FSC-H"=c(200,Inf))
 dataSubset <- filter(samp,rectGate)

}
\keyword{methods}  

