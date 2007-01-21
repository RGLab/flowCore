\name{transform}
\alias{transform}
\alias{transform,flowFrame,linearTransformation}
\alias{transform,flowFrame,quadratictransformation}
\alias{transform,flowFrame,lnTransformation}
\alias{transform,flowFrame,logTransformation}

\title{Apply transformation on a data set}
\description{transform is a function to perform
  transformations on parameters to provide suitable
visualization. 
}

\usage{
transform(`_data`,...)
}
\arguments{
  \item{`_data`}{Object of class \code{\link[flowCore]{flowFrame}}.}
   \item{...}{transformation function to be applied such as \code{\link[flowCore]{linearTransformation}}
}

\details{
 Note that the arguments on which to apply the transformation should be
 called with single quote, e.g.,transform(`_data`,`FLH-1`=2*`FLH-1`).
}

\value{Object of class \code{\link[flowCore]{flowFrame}.}
}  

\seealso{\code{\link[flowCore::transformation-class]{transformation-class}}
}

\author{N. Le Meur}
\examples{
\dontrun{
  samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  linear <- linearTransformation(transformationId="Linear-transformation", a=2, b=0)
  dataTransform <- transform(samp,`FSC-H`=linear(`FSC-H`))


  lnTransform <- lnTransformation(transformationId="LN-transformation",r=1000, d=3)
  dataTransform2 <- transform(samp,`SSC-H`=lnTransform(`SSC-H`),`FSC-H`=linear(`FSC-H`))
 }
}
\keyword{methods}  

