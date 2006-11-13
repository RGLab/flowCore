\name{filterResult-class}
\docType{class}
\alias{filterResult-class}
\alias{filterResult-accessors}
\alias{filterResult}

\title{'filterResult': a class for representing the results of the
  applyFilter method}


\description{In flow-cytometry analysis, regions in n-dimensional
  projections of the data space often have to be selected. Objects of
  this class can store the results of these selections.}

\section{Creating Objects}{
Objects can be created using methods of the generic function 
\code{  new("filterResult",}\cr
\code{    subSet= ....,   # numeric vector}\cr
\code{    gateDetails=...., # object of class character}\cr
}}


\section{Slots}{
  \describe{
    \item{\code{subSet}:}{Numeric \code{vector} of 1 and 0 to indicates
      if the data points belongs or not to the selection.}
    \item{\code{gateDetails}}{Object of class \code{character} giving
      some information of the filtering procedure.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{plot:}}{plot the results}
  }
}


\author{P.D. Haaland N. Le Meur}

\keyword{classes}