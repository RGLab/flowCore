\name{parameters}
\docType{methods}
\alias{parameters}
\alias{parameters,flowFrame,missing-method}
\alias{parameters,filter-method}
\alias{parameters,parameterFilter-method}
\alias{parameters,setOperationFilter-method}
\alias{parameters,filterReference-method}
\alias{parameters,parameterTransform-method}

\title{Obtain information about parameters for flow cytometry objects.}

\description{
  
  Many different objects in \code{flowCore} are associated with one or
  more parameters. This includes \code{\linkS4class{filter}},
  \code{\linkS4class{flowFrame}} and
  \code{\linkS4class{parameterTransform}} objects that all either
  describe or use parameters.

}

\usage{
parameters(object)
}

\arguments{

  \item{object}{Object of class \code{\linkS4class{filter}},
    \code{\linkS4class{flowFrame}} or
    \code{\linkS4class{parameterTransform}}.}

}

\value{
  
  When applied to a \code{flowFrame} object, the result is an
   \code{\linkS4class{AnnotatedDataFrame} describing the parameters
   recorded by the cytometer. For other objects it will usually return a
   vector of names used by the object for its calculations.
  
}


\section{Methods}{
  
  \describe{

    \item{object = "filter"}{Returns for all objects that inherit from
    \code{\linkS4class{filter}} a vector of parameters on which a gate
    is defined.}

    \item{object = "parameterFilter"}{see above}
    
    \item{object = "setOperationFilter"}{see above}
    
    \item{object = "filterReference"}{see above}
    
    \item{object = "flowFrame"}{Returns an
      \code{\linkS4class{AnnotatedDataFrame}} containing detailed
      descriptions about the measurement parameters of the
      \code{\linkS4class{flowFrame}}.}
    
    \item{object = "parameterTransform"}{Returns a vector of parameters
      that are transformed.}

  }
}


\author{B. Ellis, N. Le Meur, F. Hahne}

\examples{
 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 parameters(samp)
 print(samp@parameters@data)
}

\keyword{methods}

