\name{quadraticTransform}
\alias{quadraticTransform}

\title{Define quadratic transformation function to apply a data set}
\description{Quadratic transformation to perform on flowFrame parameters to provide suitable
visualization. 
}

\usage{
quadraticTransform(transformationId,a,b,c)
}
\arguments{
  \item{transformationId}{character vector to annotate the
    transformation applied}
  \item{a}{numeric multiplicative value for the quadratic factor a the
  function}
  \item{b}{numeric multiplicative factor}
  \item{c}{numeric additive factor}
}

\details{
  a*x^2 + b*x + c
}

\value{fucntion}
}  

\seealso{\code{\link[flowCore::transform]{transform}}
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

