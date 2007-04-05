\name{norm2Filter-class}
\docType{class}
\alias{norm2Filter-class}
\alias{norm2Filter}

\title{Class "norm2Filter"}
\description{Define a filter that fits a bivariate normal distribution into a data set of
  paired values and selects data points according to their standard
  deviation from the fitted distribution.}

\section{Slots}{
	 \describe{
    \item{\code{method}:}{One of \code{covMcd} or \code{cov.rob} defining method
    used for computation of covariance matrix.}
    \item{\code{scale.factor}:}{Numeric vector giving factor of standard deviations
    used for data selection (all points within \code{scalefac} standard
    deviations are selected).}
    \item{\code{transformation}:}{Object of class \code{"list"} of given
      transformation, if applicable they are applied to the data before filtering}
    \item{\code{n}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{filterId}:}{Object of class \code{"character"} referencing the filter applied}
    \item{\code{parameters}:}{Object of class \code{"ANY"} describing
      the parameters used to filter the \code{flowFrame}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{filter}"}, directly.
}
\section{Methods}{
  \describe{
    \item{summarizeFilter}{\code{signature(result = "filterResult", filter = "norm2Filter")}: ... }
  }
}
\author{F. Hahne}

\seealso{
  \code{\link[MASS::cov.rob]{cov.rob}},
  \code{\link[rrcov::covMcd]{covMcd}},
  \code{\link[flowCore]{filter}}
}
\examples{
showClass("norm2Filter")
}
\keyword{classes}
