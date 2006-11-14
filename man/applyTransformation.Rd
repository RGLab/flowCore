\name{applyTransformation}
\alias{applyTransformation}
\alias{applyTransformation,fcsFrame,transformation}
\alias{applyTransformation,fcsFrame,linearTransformation}
\alias{applyTransformation,fcsFrame,quadratictransformation}
\alias{applyTransformation,fcsFrame,lnTransformation}
\alias{applyTransformation,fcsFrame,logTransformation}

\title{Apply transformation on a data set}
\description{applyTransformation is a function to perform
  transformations on parameters to provide suitable
visualization. 
}

\usage{
applyTransformation(transformation,flowSet)
}
\arguments{
  \item{transformation}{Object of class transformation}
  \item{flowSet}{Object of class \code{\link[flowCore]{fcsFrame}}.}
}

\details{
}

\value{Object of class \code{\link[flowCore]{fcsFrame}.}
}  

\seealso{\code{\link[flowCore::transformation-class]{transformation-class}}
}

\author{N. Le Meur}
\examples{

  samp <- read.FCS(system.file("data",
   "0877408774.B08", package="flowCore"))
  linear <- new("linearTransformation",transformationId="Linear-transformation",
                parameters="FSC-H", a=2, b=0)
  dataTransform <- applyTransformation(linear,samp)


  lnTransform <- new("lnTransformation",transformationId="LN-transformation",
                parameters="FSC-H", r=1000, d=3)
  dataTransform2 <- applyTransformation(lnTransform,samp)
 

}
\keyword{methods}  

