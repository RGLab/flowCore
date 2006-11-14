\name{filter-class}
\docType{class}
\alias{filter-class}
\alias{filtergate}
\alias{rectangleGate}
\alias{polygonGate}
\alias{polytopeGate}
\alias{ellipsoidGate}
\alias{norm2Filter}
\alias{decisionTreeGate}
\alias{booleanGate}
\alias{applyGate, filter-method}

\title{'filter': a class for subsetting flow-cytometry data by defining
  regions in n-dimesional projections of the data}


\description{In flow-cytometry analysis, regions in n-dimensional
  projections of the data space often have to be selected. Objects of
  this class can store the properties of these selections.}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{  new("polygonGate",}\cr
\code{    filterId= ....,   # character vector}\cr
\code{    parameters=...., # object of class character}\cr
\code{    boundaries= ....  # matrix}\cr
\code{  )}\cr
}


\section{Slots}{
  \describe{
    \item{\code{filterId}:}{Object of class \code{character} giving the name
      of the object. You can reference the object by its name for subsequent
      operations (e.g. plotting).}
    \item{\code{parameters}}{Object of class \code{character} giving the
      name of the channels used.}
    \item{\code{boundaries}}{Object of class \code{matrix} where column represents parameters
      and row coordinates}
  }
}

\section{Methods}{
  \describe{
    \item{applyFilter:}{\code{applyFilter(filter,flowSet)}\cr applies the gating of
      object \code{filter} on data objects of class \code{\link[flowCore]{fcsFrame}} 
       or \code{matrix}.}
    \item{plot:}{Plot result}}
}


\author{P.D. Haaland N. Le Meur}

\keyword{classes}