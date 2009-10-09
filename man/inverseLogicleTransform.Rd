\name{inverseLogicleTransform}
\alias{inverseLogicleTransform}
\title{
Computes the inverse of the transform defined by the 'invLogicle_transform' function}

\description{
  inverseLogicleTransform can be use to compute the inverse of the 
  Logicle transformation. The parameters w,r,d passed as inputs should match those 
 applied to transform the data using the logicleTransform function.

}
\usage{
inverseLogicleTransform(transformationId="defaultInvLogicleTransform", w = 0,
                 r = 262144, d = 5)
}


\arguments{
  \item{transformationId}{ A name to assigned to the inverse transformation. Used
  by the transform routines. } 
  \item{w}{positive number that corresponds to the width of the negative
  data range and the range of linearized data in natural log}
  \item{r}{top of the scale data value, e.g, 10000 for common 4 decade
    data or 262144 for a 18 bit data range.}
  \item{d}{breadth of the display in decades.}
}

\references{Parks D.R., Roederer M., Moore W.A.(2006)  A new "logicle" display
  method avoids deceptive effects of logarithmic scaling for low signals
  and compensated data. CytometryA, 96(6):541-51.}
\author{N. Gopalakrishnan}

\seealso{\code{\link[flowCore]{logicleTransform}} }
\examples{

data(GvHD)
samp <- read.FCS(system.file("extdata",
  "0877408774.B08", package="flowCore"))
samp <- GvHD[[1]] 
logicle  <- logicleTransform(w=2, "logicle")
transFormedData <- transform(samp, `FSC-H`=logicle(`FSC-H`))

invLogicle <- inverseLogicleTransform(w=2,"InvLogicle")
untransFormedData  <- transform(transFormedData, `FSC-H`=invLogicle(`FSC-H`))
all.equal(exprs(samp)[,"FSC-H"], exprs(untransFormedData)[,"FSC-H"])

}



  \keyword{methods}
