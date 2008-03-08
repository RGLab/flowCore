\name{coerce}
\alias{coerce}
\alias{coerce,filter,logical-method}
\alias{coerce,filterResult,logical-method}
\alias{coerce,subsetFilter,logical-method}
\alias{coerce,unionFilter,logical-method}
\alias{coerce,complementFilter,logical-method}
\alias{coerce,factor,filterResult-method}
\alias{coerce,logical,filterResult-method}
\alias{coerce,numeric,filterResult-method}
\alias{coerce,logicalFilterResult,logical-method}
\alias{coerce,randomFilterResult,logical-method}
\alias{coerce,environment,flowSet-method}
\alias{coerce,list,flowSet-method}
\alias{coerce,intersectFilter,logical-method}
\alias{coerce,filter,call-method}
\alias{coerce,subsetFilter,call-method}
\alias{coerce,intersectFiler,call-method}
\alias{coerce,unionFilter,call-method}
\alias{coerce,complementFilter,call-method}
\alias{coerce,filterReference,concreteFilter-method}
\alias{coerce,filterReference,call-method}
\alias{coerce,formula,filter-method}
\alias{coerce,character,filter-method}
\alias{coerce,name,filter-method}
\alias{coerce,call,filter-method}
\alias{coerce,filterSet,list-method}
\alias{coerce,list,filterSet-method}
\alias{coerce,flowSet,list-method}
\alias{coerce,flowSet,flowFrame-method}
\alias{coerce,flowFrame,flowSet-method}
\alias{coerce,flowFrame,filterSet-method}

\title{Force an Object to Belong to a Class}
\description{ These functions manage the relations that allow coercing an object
     to a given class.
}

\arguments{
  \item{from, to}{The classes between which 'def' performs coercion.
    (In the case of the 'coerce' function these are objects from
    the classes, not the names of the classes, but you're not
    expected to call 'coerce' directly.)}
 }

\details{
  The function supplied as the third argument is to be
  called to implement 'as(x, to)' when 'x' has class 'from'.
  Need we add that the function should return a suitable object
  with class 'to'.
}


\author{B. Ellis}
\examples{
 samp1 <- read.FCS(system.file("extdata","0877408774.E07", package="flowCore"))
 samp2 <- read.FCS(system.file("extdata","0877408774.B08",package="flowCore"))
 samples <-list("sample1"=samp1,"sample2"=samp2)
 experiment <- as(samples,"flowSet")
}
\keyword{methods}  

