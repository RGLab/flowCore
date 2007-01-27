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
transform(`_data`,\dots)
}
\arguments{
  \item{_data}{Object of class \code{\link[flowCore]{flowFrame}}.}
   \item{\dots}{transformation function to be applied such as \code{\link[flowCore]{linearTransformation}}
   }
 }

 \details{
  The \code{\dots} arguments to \code{transform} function are tagged
  vector expressions, which are evaluated in the data frame
  \code{_data}.  The tags are matched against \code{names(_data)}, and for
  those that match, the value replace the corresponding variable in
  \code{_data}, and the others are appended to \code{_data}.
  
  Note that the arguments on which to apply the transformation should be
  called with single quote.
}

\value{Object of class \code{\link[flowCore]{flowFrame} with modified value.}
}  

\seealso{\code{\link[flowCore::transformation-class]{transformation-class}}
}

\author{N. Le Meur}
\examples{
\dontrun{
  samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))
  linear <- linearTransform(transformationId="Linear-transformation", a=2, b=0)
  dataTransform <- transform(samp,`FSC-H`=linear(`FSC-H`))


  lnTransform <- lnTransform(transformationId="LN-transformation",r=1000, d=3)
  dataTransform2 <- transform(samp,`SSC-H`=lnTransform(`SSC-H`),`FSC-H`=linear(`FSC-H`))
 }
}
\keyword{methods}  

