\name{filterSet-class}
\docType{class}
\alias{filterSet-class}
\alias{filterSet}
\alias{filterReference,filterSet,character-method}
\alias{identifier,filterSet-method}
\alias{identifier,filterSet,character-method}
\alias{names,filterSet-method}
\alias{show,filterSet-method}
\alias{sort,filterSet-method}
\alias{\%subset\%,filterSet,filter-method}
\alias{[,filterSet,character-method}
\alias{[[,filterSet,character-method}
\alias{[[<-,filterSet,NULL,ANY,filter-method}
\alias{[[<-,filterSet,NULL,ANY,formula-method}
\alias{[[<-,filterSet,character,ANY,formula-method}
\alias{[[<-,filterSet,character,ANY,filter-method}
\alias{[[<-,filterSet,missing,ANY,filter}
\alias{[[<-,filterSet,ANY,ANY,filterReference-method}
\alias{[[<-,filterSet,missing,ANY,filter-method}

\title{Class filterSet}

\description{A container for a collection of related filters. }

\section{Objects from the Class}{
  
  There are several ways to create a \code{filterSet} object. There is
  the \code{\link{filterSet}} constructor, which creates an empty
  \code{filterSet} object (see the details section for more
  information). \code{filterSet} objects can also be coerced to and from
  \code{list} objects using the \code{as} function.

}

\section{Slots}{
  \describe{
    \item{\code{env}:}{The environment that actually holds the filters}
    
    \item{\code{name}:}{A more descriptive name of the set.}
    
  }
}

\section{Methods}{
 \describe{
  \item{names}{An unsorted list of the names of the filters contained
    within the set.}
 
  \item{sort}{Returns a topological sort of the names of the filters contained within the 
set. Primarily used by internal functions (such as \code{\link{filter}}), this method is 
also useful for planning gating strategy layouts and the like.}

  \item{filterReference}{Retrieves references to a filter inside a
    filterSet}
  
  \item{[}{Returns the filter reference used inside the filter. See
    Details.}
  
  \item{[[}{Retrieves the actual filters from a filterSet. Note that
    composed filters can still contain references.}
  
  \item{[[<-}{Put a filter into a filterSet. As a convenience, assigning
    to the ``""'' or \code{NULL} name will use the filter's name for
    assignment. Composed filters can be added easily using formulas
    rather than attempting to construct filters the long way. The
    formula interface is also lazy, allowing you to add filters in any
    order.}
}
}

\section{Details}{
  
  \code{filterSet} objects are intended to provide a convenient grouping
  mechanism for a particular gating strategy. To accomplish this, much
  like the \code{flowSet} object, the \code{filterSet} object introduces
  reference semantics through the use of an \code{environment}, allowing
  users to change an upstream filter via the usual assignment mechanism
  and have that change reflected in all dependent filters. We do this by
  actually creating two filters for each filter in the \code{filterSet}.
  The first is the actual concrete filter, which is assigned to a
  variable of the form \code{.name} where \code{name} is the original
  filter name. A second \code{filterReference} filter is the created
  with the original name to point to the internal name. The allows us to
  evaluate a \code{formula} in the environment without creating a copy
  of the original filter.

}

\author{B. Ellis}

\seealso{
  
  \code{\link{filterSet}}

}

\examples{

fs = new("filterSet")

## Simple assignment. Note that the filterId slot for the rectangle gate
## is changed.
fs[["filter1"]] = rectangleGate("FSC-H"=c(.2,.8),"SSC-H"=c(0,.8))

## Convenience assignment using the filterId slot.
fs[[""]] = norm2Filter("FSC-H","SSC-H",scale.factor=2,filterId="Live Cells")

## We also support formula interfaces. These two statements are equivalent.
fs[["Combined"]] = ~ filter1 \%subset\% `Live Cells`
fs[[""]] = Combined ~ filter1 \%subset\% `Live Cells`
fs
as(fs,"list")
as(as(fs,"list"),"filterSet")

}
\keyword{classes}
