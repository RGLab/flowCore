\name{logicleTransform}
\alias{logicleTransform}

\title{Computes a transform using the 'logicle_transform' function}
\description{
  
  Logicle transformation creates a subset of \code{\link{biexponentialTransform}}
  hyperbolic sine transformation functions that provides several advantages over linear/log 
  transformations for display of flow cytometry data.
}
\usage{
logicleTransform(transformationId="defaultLogicleTransform", w = 0,
                 t = 262144, m = 4.5, a=0, tol =.Machine$double.eps^0.8,
                 maxit = as.integer(5000)) 
}

\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used
  by the transform/filter  routines. } 
  \item{w}{ w is the linearization width in asymptotic decades. W should be >= 0 and 
  determines the slope of  transformation at zero.
  w can be estimated using the equation 
  w=m-log10(t/abs(r)))/2 , where r is the most negative value to be included in the 
  display
  }
  \item{t}{Top of the scale data value, e.g, 10000 for common 4 decade
    data or 262144 for a 18 bit data range. t should be greater than zero}
  \item{m}{m is the full width of the transformed display in asymptotic decades. m
  should be greater than zero}
  \item{a}{Additional negative range to be included in the display in asymptotic 
   decades. Positive values of the argument brings additional negative input values
   into the transformed display viewing area. Default value is zero corresponding to a Standard logicle function.}
  \item{tol}{Acceptable tolerance of the root value}
  \item{maxit}{ Maximum iterations allowed for the root search process}
}

\references{Parks D.R., Roederer M., Moore W.A.(2006)  A new "logicle" display
  method avoids deceptive effects of logarithmic scaling for low signals
  and compensated data. CytometryA, 96(6):541-51.}
\author{ B. Ellis N. LeMeur, N Gopalakrishnan}

\seealso{\code{\link[flowCore]{biexponentialTransform}} }
\examples{
data(GvHD)
samp <- read.FCS(system.file("extdata",
  "0877408774.B08", package="flowCore"))
samp <- GvHD[[1]] 
logicle  <- logicleTransform(w=1, "logicle")
after <- transform(samp, `FSC-H`=logicle(`FSC-H`))
  }

  \keyword{methods}
