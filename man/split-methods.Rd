\name{split-methods}
\docType{methods}
\alias{split-methods}
\alias{split}
\alias{split,flowFrame,ANY-method}
\alias{split,flowFrame,factor-method}
\alias{split,flowFrame,character-method}
\alias{split,flowFrame,numeric-method}
\alias{split,flowFrame,filter-method}
\alias{split,flowFrame,filterSet-method}
\alias{split,flowFrame,logicalFilterResult-method}
\alias{split,flowFrame,manyFilterResult-method}
\alias{split,flowFrame,multipleFilterResult-method}
\alias{split,flowSet,ANY-method}
\alias{split,flowSet,character-method}
\alias{split,flowSet,factor-method}
\alias{split,flowSet,list-method}
\alias{split,flowSet,numeric-method}
\alias{split,flowSet,filter-method}


\title{ Methods to split flowFrames and flowSets according to filters }

\description{
  
  Divide a flow cytometry data set into several subset according to the
  results of a filtering operation.

}

\section{Methods}{

  \code{\link{flowFrame}} methods:

  \describe{
  
    \item{x = "flowFrame", f = "ANY", drop = "ANY"}{ Catch all input and
      cast an error if there is no method for \code{f} to dispatch to. }

    \item{x = "flowFrame", f = "factor", drop = "ANY"}{ Split a
      \code{\link{flowFrame}} by a factor variable. Length of \code{f}
      needs to the same as \code{nrow(x)}. }

    \item{x = "flowFrame", f = "character", drop = "ANY"}{ Coerce
      \code{f} to a factor and split on that. }

    \item{x = "flowFrame", f = "numeric", drop = "ANY"}{ Coerce \code{f} to
      a factor and split on that. }

    \item{x = "flowFrame", f = "filter", drop = "ANY"}{ First applies
      the \code{\linkS4class{filter}} to the
      \code{\linkS4class{flowFrame}} and then splits on the resulting
      \code{\linkS4class{filterResult}} object. }

    \item{x = "flowFrame", f = "filterSet", drop = "ANY"}{ First applies
      the \code{\linkS4class{filterSet}} to the
      \code{\linkS4class{flowFrame}} and then splits on the resulting
      final \code{\linkS4class{filterResult}} object. }

    \item{x = "flowFrame", f = "logicalFilterResult", drop = "ANY"}{ Split
      into the two subpopulations (in and out of the gate). }

    \item{x = "flowFrame", f = "manyFilterResult", drop = "ANY"}{ Split
      into the several subpopulations identified by the filtering
      operation. Instead of returning a list, the additional logical
      argument code{flowSet} makes the method return an object od class
      \code{\linkS4class{flowSet}}. The \code{population} argument takes
      a character vector indicating the subpopulations to use for
      splitting (as identified by the population name in the
      \code{filterDetails} slot).}

    \item{x = "flowFrame", f = "multipleFilterResult", drop = "ANY"}{
      Split into the several subpopulations identified by the filtering
      operation. Instead of returning a list, the additional logical
      argument code{flowSet} makes the method return an object od class
      \code{\linkS4class{flowSet}}. The \code{population} argument takes
      a character vector indicating the subpopulations to use for
      splitting (as identified by the population name in the
      \code{filterDetails} slot). }

    }

    \code{\linkS4class{flowSet}} methods:

     \describe{

    \item{x = "flowSet", f = "ANY", drop = "ANY"}{ Catch all input and cast
      an error if there is no method for \code{f} to dispatch to.  }

    \item{x = "flowSet", f = "factor", drop = "ANY"}{ Split a
      \code{\link{flowSet}} by a factor variable. Length of \code{f}
      needs to the same as \code{length(x)}. }

    \item{x = "flowSet", f = "character", drop = "ANY"}{ Coerce \code{f} to
      a factor and split on that. }

    \item{x = "flowSet", f = "numeric", drop = "ANY"}{ Coerce \code{f}
      to a factor and split on that. }

    \item{x = "flowSet", f = "list", drop = "ANY"}{ Split a
      \code{\link{flowSet}} by a list of
      \code{\linkS4class{filterResult}}s (as typically returned by
      filtering operations on a \code{\linkS4class{flowSet}}). The
      length of the list has to be equal to the length of the
      \code{\linkS4class{flowSet}} and every list item needs to be a
      \code{\linkS4class{filterResult}} of equal class.  Instead of
      returning a list, the additional logical argument code{flowSet}
      makes the method return an object od class
      \code{\linkS4class{flowSet}}. The \code{population} argument takes
      a character vector indicating the subpopulations to use for
      splitting (as identified by the population name in the
      \code{filterDetails} slot). }
    
  }
}

\details{The flow data structures in \code{flowCore} can be split into
  subsets on various levels:

  \code{\linkS4class{flowFrame}}: row-wise splitting of the raw data. In
  most cases, this will be done according to the outcome of a filtering
  operation, either using a filter that idenitfiers more than one
  sub-poupulation or by a logical filter, in which case the data is
  split into "in the filter" and "not in the filter". In addtion, the
  data can be split according to a factor (or a numeric or character
  vector that can be coerced into a factor).

  code{\linkS4class{flowSet}}: can be either split into subsets of
  \code{\linkS4class{flowFrame}}s according to a factor or a vector that
  can be coerced into a factor, or each individual
  \code{\linkS4class{flowFrame}} into subpopulations based on the
  \code{\linkS4class{filter}}s or \code{\linkS4class{filterResult}}s
  provided as a list of equal length.

  In both cases, if the object to split inherits from class
  \code{\linkS4class{multipleFilterResult}} or
  \code{\linkS4class{manyFilterResult}}, the optional argument
  \code{population} can be provided to only select certain
  subpopulations. These subpopulations are identified by their
  \code{population} name stored in the \code{filterDetails} slot.

  The default return value for any of the \code{split} methods is a
  list, but the optional logical argument \code{filterSet} can be used
  to return a \code{\linkS4class{flowSet}} instead.
  
}

\author{F Hahne, B. Ellis, N. Le Meur}

\examples{

## split a flowFrame
data(GvHD)
samp <- GvHD[[1]]
qGate <- quadGate(filterId="qg", "FSC-H"=200, "SSC-H"=400)
fres <- filter(samp, qGate)
split(samp, qGate)
split(samp, fres, flowSet=TRUE)

## only keep one population
split(samp, fres, population="tl")


## split the whole set
split(GvHD, qGate, flowSet=TRUE, population=c("tl", "br"))

## now according to a factor
split(GvHD, pData(GvHD)$Patient)

}



\keyword{methods}
