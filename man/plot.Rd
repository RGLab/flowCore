\name{plot}
\alias{plot}
\alias{plot,flowFrame}

\title{Plot floFrame object}
\description{
  A density plot of raw flowcytometry data or a scatterplot of
  filterResult object, i.e. data of filtered the flowFrame object.
}
\usage{
plot(x, y, filter=NULL,...)
}
\arguments{
  \item{x}{flowFrame object}
  \item{y}{missing or character vector}
  \item{...}{\code{plot} options}
}

\author{P. Haaland and N. Le Meur}

\examples{
samp <- read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
plot(samp)

}
\keyword{methods}