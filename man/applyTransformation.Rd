\name{applyTransformation}
\alias{applyTransformation}
\alias{applyTransformation,flowFrame,transformation}
\alias{applyTransformation,flowFrame,linearTransformation}
\alias{applyTransformation,flowFrame,quadratictransformation}
\alias{applyTransformation,flowFrame,lnTransformation}
\alias{applyTransformation,flowFrame,logTransformation}

\title{Apply transformation on a data set}
\description{applyTransformation is a function to perform
  transformations on parameters to provide suitable
visualization. 
}

\usage{
applyTransformation(transformation,flowObject)
}
\arguments{
  \item{transformation}{Object of class transformation}
  \item{flowObject}{Object of class \code{\link[flowCore]{flowFrame}}.}
}

\details{
}

\value{Object of class \code{\link[flowCore]{flowFrame}.}
}  

\seealso{\code{\link[flowCore::transformation-class]{transformation-class}}
}

\author{N. Le Meur}
\examples{

  samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  linear <- new("linearTransformation",transformationId="Linear-transformation",
                parameters="FSC-H", a=2, b=0)
  dataTransform <- applyTransformation(linear,samp)


  lnTransform <- new("lnTransformation",transformationId="LN-transformation",
                parameters="FSC-H", r=1000, d=3)
  dataTransform2 <- applyTransformation(lnTransform,samp)
 

}
\keyword{methods}  

