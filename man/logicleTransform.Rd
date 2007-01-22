\name{logicleTransformation}
\alias{logicleTransformation}


\title{Define logicle transformation function to apply a data set}
\description{Logicle transformation to perform on flowFrame parameters to provide suitable
visualization. 
}

\usage{
logicleTransform(w=0,r=262144,d=5,...)
}
\arguments{
  \item{transformationId}{character string to annotate the
    transformation applied}
  \item{r}{positive numeric value corresponding to the analog-to-digital resolution}
  \item{d}{positive numeric value corresponding to the number of decades for the dynamic
range of parameter.}
  \item{w}{positive numeric value}
  
}

\details{
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
  logicle <- logicleTransform(transformationId="Logicle-transformation", w=0,r=262144,d=5)
  dataTransform <- transform(samp,`FSC-H`=logicle(`FSC-H`))

 }
}
\keyword{methods}  

