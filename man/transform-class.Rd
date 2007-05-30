\name{transform-class}
\docType{class}
\alias{transform}
\alias{transform,missing-method}
\alias{transform-class}
\alias{summary,transform-method}

\title{'transform': a class for transforming flow-cytometry data by
  applying scale factors.}
\description{
	Transform objects are simply functions that have been extended to allow
for specialized dispatch. All of the ``...Transform'' constructors return 
functions of this type for use in one of the transformation modalities.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"function"}}
  }
}
\section{Methods}{
  \item{\code{summary}}{Return the parameters}
}

\author{N LeMeur}

\seealso{\code{\link[flowCore]{linearTransform}},
  \code{\link[flowCore]{lnTransform}},
  \code{\link[flowCore]{logicleTransform}},
  \code{\link[flowCore]{biexponentialTransform}},
  \code{\link[flowCore]{arcsinhTransform}},
  \code{\link[flowCore]{quadraticTransform}},  \code{\link[flowCore]{logTransform}}}

\examples{

cosTransform <- function(transformId, a=1, b=1){
  t = new("transform", .Data = function(x) cos(a*x+b))
  t@transformationId = transformId
  t
}

cosT <- cosTransform(transformId="CosT",a=2,b=1)

summary(cosT)
}
\keyword{classes}
