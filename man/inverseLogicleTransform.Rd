\name{inverseLogicleTransform}
\alias{inverseLogicleTransform}
\title{
Computes the inverse of the transform defined by the 'invLogicle_transform' function}

\description{
  inverseLogicleTransform can be use to compute the inverse of the 
  Logicle transformation. The parameters w,t,m,a passed as inputs should match those 
 applied to transform the data using the logicleTransform function.

}
\usage{
inverseLogicleTransform(transformationId="defaultInvLogicleTransform", w = 0,
                 t = 262144, m = 4.5,a=0)
}


\arguments{
  \item{transformationId}{ A name to assigned to the inverse transformation. Used
  by the transform routines. } 
  \item{w}{w is the linearization width in asymptotic decades. W should be >= 0 and 
  determines the slope of  transformation at zero}
  \item{t}{Top of the scale data value, e.g, 10000 for common 4 decade
    data or 262144 for a 18 bit data range. t should be greater than zero.}
  \item{m}{m is the full width of the transformed display in asymptotic decades. m
  should be greater than zero.}
   \item{a}{Additional negative range to be included in the display in asymptotic 
   decades. Positive values of the argument brings additional negative input values
   into the transformed display viewing area. Default value is zero corresponding to a Standard logicle function.}

}

\references{Parks D.R., Roederer M., Moore W.A.(2006)  A new "logicle" display
  method avoids deceptive effects of logarithmic scaling for low signals
  and compensated data. CytometryA, 96(6):541-51.}
\author{N. Gopalakrishnan}

\seealso{\code{\link[flowCore]{logicleTransform}} }
\examples{

data(GvHD)
samp <- GvHD[[1]] 
logicle  <- logicleTransform(w=2, "logicle")
transFormedData <- transform(samp, `FSC-H`=logicle(`FSC-H`))

invLogicle <- inverseLogicleTransform(w=2,"InvLogicle")
untransFormedData  <- transform(transFormedData, `FSC-H`=invLogicle(`FSC-H`))
all.equal(exprs(samp)[,"FSC-H"], exprs(untransFormedData)[,"FSC-H"])

}
\keyword{methods}
