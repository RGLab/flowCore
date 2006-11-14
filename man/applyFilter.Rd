\name{applyFilter}
\alias{applyFilter}
\alias{applyFilter,fcsFrame,filter}
\alias{applyFilter,fcsFrame,rectangleGate}
\alias{applyFilter,fcsFrame,polygonGate}
\alias{applyFilter,fcsFrame,norm2Filter}

\title{Apply filter on a data set}
\description{applyFilter is a function for subsetting data by applying}

\usage{
applyFilter(filter,flowSet,parent)
}
\arguments{
  \item{filter}{Object of class gate and slots parameters and boundaries}
  \item{flowSet}{Object of class \code{\link[flowCore]{fcsFrame}}.}
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
 samp <- read.FCS(system.file("data","0877408774.B08", package="flowCore"))
 limitMin <- 200
 limitMax<-  Inf
 filter <- new("rectangleGate",filterId="nonDebris",parameters="FSC-H",min= limitMin,max=limitMax)
 dataSubset <- applyFilter(filter,samp)

}
\keyword{methods}  

