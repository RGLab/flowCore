\name{linearTransform}
\alias{linearTransform}

\title{Create the definition of linear transformation function to be applied on a data set }
\description{Create the definition of linear Transformation that will be
  applied on some parameter via the \code{transform} method. 
  The definition of this function is currently  x <- a*x+b
}
\usage{
linearTransform(transformationId, a = 1, b = 0)
}

\arguments{
  \item{transformationId}{character stinger to identifier ~~Describe \code{transformationId} here~~ }
  \item{a}{double that correponds to the multiplicative factor in the equation}
  \item{b}{double that correponds to the additive factor in the equation}
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\author{ N. Le Meur }

\seealso{ \code{\link{transform-class}}, \code{\link{transform}}}
\examples{
\dontrun{samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  linear <- linearTransform(transformationId="Linear-transformation", a=2, b=0)
  dataTransform <- transform(samp,`FSC-H`=linear(`FSC-H`))
}

}

\keyword{methods}

