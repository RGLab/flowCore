\name{transformation-class}
\docType{class}
\alias{transformation-class}
\alias{linearTransformation}
\alias{quadraticTransformation}
\alias{lnTransformation}
\alias{logTransformation}
\alias{applyTransformation, transformation-method}

\title{'Transformation': a class for transforming flow-cytometry data by
  applying scale factors.}


\description{While analyzing flow cytometry data, various parameter transformations are performed to provide suitable visualization.}

\section{Creating Objects}{
Objects can be created using methods of the generic function. Here is an
example for a linear transformation of the type newFCS =a*parameters + b 
\code{  new("linearTransformation",}\cr
\code{    transformationId= ....,   # character vector}\cr
\code{    parameters=...., # object of class character}\cr
\code{    a= ....  # factor}\cr
\code{    b= ....  # factor}\cr
\code{  )}\cr
}


\section{Slots}{
  \describe{
    \item{\code{transformationId}:}{Object of class \code{character} giving the name
      of the object. You can reference the object by its name for subsequent
      operations (e.g. plotting).}
    \item{\code{parameters}}{Object of class \code{character} or \code{numeric}
      giving the name of the channels used.}
    \item{\code{a}}{Object of class \code{numeric}, scaling factor.}
    \item{\code{b}}{Object of class \code{numeric}, scaling factor.}
  }
}

\section{Methods}{
  \describe{
    \item{applyTransformation:}{\code{applyTransformation(transformation,flowObject)}\cr applies the transformation on data objects of class \code{\link[flowCore]{flowFrame}}.}
    \item{plot:}{Plot result}
  }
}


\author{N. Le Meur}

\keyword{classes}