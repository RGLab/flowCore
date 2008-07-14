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
  results of a filtering operation. There are also methods available to
  split according to a factor variable.

}

\section{Methods}{

  \code{\link{flowFrame}} methods:

  \describe{
  
    \item{x = "flowFrame", f = "ANY", drop = "ANY"}{ Catch all input and
      cast an error if there is no method for \code{f} to dispatch to. }

    \item{x = "flowFrame", f = "factor", drop = "ANY"}{ Split a
      \code{\link{flowFrame}} by a factor variable. Length of \code{f}
      should be the same as \code{nrow(x)}, otherwise it will be
      recyled, possibly leading to undesired outcomes. The optional
      argument \code{drop} works in the usual way, in that it removes
      empty levels from the factor before splitting.}

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
      into the two subpopulations (in and out of the gate). The optional
      argument \code{population} can be used to control the names of the
      results. }

    \item{x = "flowFrame", f = "manyFilterResult", drop = "ANY"}{ Split
      into the several subpopulations identified by the filtering
      operation. Instead of returning a list, the additional logical
      argument code{flowSet} makes the method return an object of class
      \code{\linkS4class{flowSet}}. The optional \code{population}
      argument takes a character vector indicating the subpopulations to
      use for splitting (as identified by the population name in the
      \code{filterDetails} slot).}

    \item{x = "flowFrame", f = "multipleFilterResult", drop = "ANY"}{
      Split into the several subpopulations identified by the filtering
      operation. Instead of returning a list, the additional logical
      argument code{flowSet} makes the method return an object od class
      \code{\linkS4class{flowSet}}. The optional \code{population}
      argument takes a character vector indicating the subpopulations to
      use for splitting (as identified by the population name in the
      \code{filterDetails} slot). }

    }

    \code{\linkS4class{flowSet}} methods:

     \describe{

    \item{x = "flowSet", f = "ANY", drop = "ANY"}{ Catch all input and cast
      an error if there is no method for \code{f} to dispatch to.  }

    \item{x = "flowSet", f = "factor", drop = "ANY"}{ Split a
      \code{\link{flowSet}} by a factor variable. Length of \code{f}
      needs to be the same as \code{length(x)}. The optional argument
      \code{drop} works in the usual way, in that it removes empty
      levels from the factor before splitting. }

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
      \code{\linkS4class{filterResult}} of equal class with the same
      parameters.  Instead of returning a list, the additional logical
      argument code{flowSet} makes the method return an object of class
      \code{\linkS4class{flowSet}}. The optional \code{population}
      argument takes a character vector indicating the subpopulations to
      use for splitting (as identified by the population name in the
      \code{filterDetails} slot). Note that this also implies common
      population names for all\code{\linkS4class{filterResult}}s in the
      list. }
    
  }
}

\details{

  The splitting operation in the context of
  code{\linkS4class{flowFrame}}s and code{\linkS4class{flowSet}}s is the
  logical extension of subsetting. While the latter only returns the
  events contained within a gate, the former splits the data into the
  groups of events cotained within and those not contained within a
  particular gate. This concept is extremely useful in applications
  where gates describe the distinction between positivity and negativity
  for a particual marker.


  The flow data structures in \code{flowCore} can be split into
  subsets on various levels:

  \code{\linkS4class{flowFrame}}: row-wise splitting of the raw data. In
  most cases, this will be done according to the outcome of a filtering
  operation, either using a filter that idenitfiers more than one
  sub-poupulation or by a logical filter, in which case the data is
  split into two populations: "in the filter" and "not in the
  filter". In addtion, the data can be split according to a factor (or a
  numeric or character vector that can be coerced into a factor).

  \code{\linkS4class{flowSet}}: can be either split into subsets of
  \code{\linkS4class{flowFrame}}s according to a factor or a vector that
  can be coerced into a factor, or each individual
  \code{\linkS4class{flowFrame}} into subpopulations based on the
  \code{\linkS4class{filter}}s or \code{\linkS4class{filterResult}}s
  provided as a list of equal length.

  Splitting has a special meaning for filters that result in
  \code{\linkS4class{multipleFilterResult}}s or
  \code{\linkS4class{manyFilterResult}}s, in which case simple
  subsetting doesn't make much sense (there are multiple populations
  that are defined by the gate and it is not clear which of those should
  be used for the subsetting operation). Accordingly, splitting of
  multipleFilterResults creates multiple subsets. The argument
  \code{population} can be used to limit the ouput to only one or some
  of the resulting subsets. It takes as values a character vector of
  names of the populations of interest. See the documentation of the
  different filter classes on how population names can be defined and
  the respective default values. For splitting of
  \code{\linkS4class{logicalFilterResult}}s, the \code{population}
  argument can be used to set the population names since there is no
  reasonable default other than the name of the gate. The content of the
  argument \code{prefix} will be prepended to the population names and
  '+' or '-' are finally appended allowing for more flexible naming
  schemes.


  The default return value for any of the \code{split} methods is a
  list, but the optional logical argument \code{filterSet} can be used
  to return a \code{\linkS4class{flowSet}} instead. This only applies
  when splitting \code{\linkS4class{flowFrame}}s, splitting of
  \code{\linkS4class{flowSet}}s always results in lists of
  \code{\linkS4class{flowSet}} objects.
  
}

\author{F Hahne, B. Ellis, N. Le Meur}

\examples{

data(GvHD)
qGate <- quadGate(filterId="qg", "FSC-H"=200, "SSC-H"=400)

## split a flowFrame by a filter that creates
## a multipleFilterResult
samp <- GvHD[[1]]
fres <- filter(samp, qGate)
split(samp, qGate)

## return a flowSet rather than a list
split(samp, fres, flowSet=TRUE)

## only keep one population
names(fres)
split(samp, fres, population="FSC-Height+SSC-Height+")


## split the whole set, only keep two populations
split(GvHD, qGate, population=c("FSC-Height+SSC-Height+",
"FSC-Height-SSC-Height+"))

## now split the flowSet according to a factor
split(GvHD, pData(GvHD)$Patient)

}



\keyword{methods}
