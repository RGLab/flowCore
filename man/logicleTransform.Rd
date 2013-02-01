\name{logicleTransform}
\alias{logicleTransform}
\alias{estimateLogicle}
\title{Computes a transform using the 'logicle_transform' function}
\description{
  
  Logicle transformation creates a subset of \code{\link{biexponentialTransform}}
  hyperbolic sine transformation functions that provides several advantages over linear/log 
  transformations for display of flow cytometry data. (The logicleTransform method 
  makes use of the C++ implementation of the logicle transform contributed
  by Wayne Moore et al.) 
}
\usage{
logicleTransform(transformationId="defaultLogicleTransform", w = 0.5, t = 262144,
        m = 4.5, a = 0)
 estimateLogicle(x, channels,...) 
}

\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used
  by the transform/filter  routines. } 
  \item{w}{ w is the linearization width in asymptotic decades. w should be > 0 and 
  determines the slope of  transformation at zero.
  w can be estimated using the equation 
  w=(m-log10(t/abs(r)))/2, where r is the most negative value to be included in the 
  display
  }
  \item{t}{Top of the scale data value, e.g, 10000 for common 4 decade
    data or 262144 for a 18 bit data range. t should be greater than zero}
  \item{m}{m is the full width of the transformed display in asymptotic decades. m
  should be greater than zero}
  \item{a}{Additional negative range to be included in the display in asymptotic 
   decades. Positive values of the argument brings additional negative input values
   into the transformed display viewing area. Default value is zero corresponding 
   to a Standard logicle function.}
  \item{x}{ Input flow frame for which the logicle transformations are to be estimated.}
  \item{channels}{ Flouresence channels for which the logicle transformation is to 
      be estimated.}
}

\references{Parks D.R., Roederer M., Moore W.A.(2006)  A new "logicle" display
  method avoids deceptive effects of logarithmic scaling for low signals
  and compensated data. CytometryA, 96(6):541-51.}
\author{ Wayne Moore, N Gopalakrishnan}

\seealso{\code{\link[flowCore]{inverseLogicleTransform}}, 
    \code{\link[flowCore]{estimateLogicle} }}
\examples{
data(GvHD)
samp <- GvHD[[1]] 
## User defined logicle function
lgcl <- logicleTransform( w = 0.5, t= 10000, m =4.5)
after <- transform(samp, `FL1-H` = lgcl(`FL1-H`), `FL2-H` = lgcl(`FL2-H`))
invLgcl <- inverseLogicleTransform(trans = lgcl)
before <- transform (after,  `FL1-H`= invLgcl(`FL1-H`), `FL2-H` = invLgcl(`FL2-H`))

## Automatically estimate the logicle transformation based on the data
lgcl <- estimateLogicle(samp, channels = c("FL1-H", "FL2-H", "FL3-H", "FL2-A", "FL4-H"))
## transform  parameters using the estimated logicle transformation
after <- transform(samp, lgcl)

}

\keyword{methods}
