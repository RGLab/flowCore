\name{linearTransform}
\alias{linearTransform}


\title{Define a linear transformation function to apply a data set}
\description{Linear transform to perform on flowFrame parameters to provide suitable
visualization. 
}

\usage{
linearTransform(transformationId,a,b)
}
\arguments{
  \item{transformationId}{character vector to annotate the
    transformation applied}
  \item{a}{numeric multiplicative factor}
  \item{b}{numeric additive factor}
}

\details{
  a*x+b
}

\value{function}
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

 }
}
\keyword{methods}  

