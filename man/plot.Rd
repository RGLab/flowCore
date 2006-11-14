\name{plot}
\alias{plot}
\alias{basicplot,fcsFrame}

\title{Plot fcsframe object}
\description{
  A density plot of raw flowcytometry data or a scatterplot of
  filterResult object, i.e. data of filtered the fcsframe object.
}
\usage{
basicplot(x, y, filter=NULL, plotParameters=c(1,2),parent,colParent="Grey",colSubSet="Blue",xlim,ylim,...)
}
\arguments{
  \item{x}{fcsframe object}
  \item{y}{filterResult Object}
  \item{filter}{gate Object}
  \item{plotParameters}{parameters to plot}
  \item{parent}{}
  \item{colParent}{}
  \item{colSubSet}{}
  \item{xlim}{numerical vector of the range of the x variable
    (horizontal axis)}
  \item{ylim}{numerical vector of the range of the y variable
    (vertical axis)}
  \item{...}{\code{plot} options}
}

\author{P. Haaland and N. Le Meur}

\examples{
 samp <- read.FCS(system.file("data","0877408774.B08",package="flowCore"))
 plot(samp)

##After filtering by a triangle gate
triangleCut <- matrix(c(300,400,600,300,400,300),ncol=2,nrow=3)
triangleF <- new("polygonGate",filterId="nonDebris",
            parameters=c("FSC-H","SSC-H"),boundaries= triangleCut)
sampTrg <- applyFilter(triangleF,samp)
plot(samp,sampTrg,filter=triangleF,pch="*")
}
\keyword{methods}