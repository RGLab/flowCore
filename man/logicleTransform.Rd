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
                 r = 262144, d = 5,tol =.Machine$double.eps^0.8,
                 maxit = as.integer(5000)) 
}


\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used
  by the transform/filter integration routines. } 
  \item{w}{positive number that corresponds to the width of the negative
  data range and the range of linearized data in natural log}
  \item{r}{top of the scale data value, e.g, 10000 for common 4 decade
    data or 262144 for a 18 bit data range.}
  \item{d}{breath of the display in natural logarithm units.}
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
logicle  <- logicleTransform(w=2, "logicle")
after <- transform(samp, `FSC-H`=logicle(`FSC-H`))
  }

  \keyword{methods}
