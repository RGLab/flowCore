\name{inverseLogicleTransform}
\alias{inverseLogicleTransform}
\title{
Computes the inverse of the transform defined by the 'logicleTransform' function}

\description{
  inverseLogicleTransform can be use to compute the inverse of the 
  Logicle transformation. The parameters w, t, m, a for calculating 
  the inverse are obtained from the 'trans' input passed to the 
  'inverseLogicleTransform' function. (The inverseLogicleTransform 
  method  makes use of the C++ implementation of the inverse logicle transform 
  contributed by Wayne Moore et al.) 

}
\usage{
inverseLogicleTransform(transformationId, trans)
}


\arguments{
  \item{transformationId}{ A name to assigned to the inverse transformation. Used
  by the transform routines. } 
  \item{trans}{An object of class 'transform' created using the 'logicleTransform'
   function. The parameters w, t, m, a for calculating the inverse are obtained from the 'trans' input passed to the 
  'inverseLogicleTransform' function. }
}

\references{Parks D.R., Roederer M., Moore W.A.(2006)  A new "logicle" display
  method avoids deceptive effects of logarithmic scaling for low signals
  and compensated data. CytometryA, 96(6):541-51.}
\author{Wayne Moore, N. Gopalakrishnan}

\seealso{\code{\link[flowCore]{logicleTransform}}}
\examples{

data(GvHD)
samp <- GvHD[[1]] 
logicle  <- logicleTransform(t = 10000, w = 0.5, m = 4.5 , a =0 ,"logicle")
## transform FL1-H parameter using logicle transformation
after <- transform(samp, transformList('FL1-H', logicle))

## Inverse transform the logicle transformed data to retrieve the original data
invLogicle <- inverseLogicleTransform(trans = logicle)
before <- transform (after, transformList('FL1-H', invLogicle))

}
\keyword{methods}
