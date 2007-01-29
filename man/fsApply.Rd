\name{fsApply}
\alias{fsApply}
\alias{fsApply,flowSet,ANY}


\title{Apply a Function over values in a flowSet}

\description{'fsApply' applies 'FUN' to the named values from an
  flowSet environment and returns the results as a list. The user can request that all named
     objects are used (normally names that begin with a dot are not).
     The output is not sorted and no parent environments are searched.
}

\usage{
    fsApply(x, FUN, simplify=TRUE, use.exprs=FALSE,\dots )
}

\arguments{
  \item{x}{\code{\link[flowCore]{flowSet}} to be used}
  \item{FUN}{the function to be applied to each element of \code{x}}
  \item{simplify}{logical; should the result be simplified to a vector
    or matrix if possible?}
  \item{use.exprs}{logical: should the \code{FUN} be applied on the
    \code{\link[flowCore]{flowFrame}} object}
  \item{\dots}{optional arguments to \code{FUN}.}
  
}
\seealso{
  \code{\link{apply}}
}
\examples{
fcs.loc <- system.file("extdata",package="flowCore")
file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")

samp <- read.flowSet(file.location[1:3]) 
fsApply(samp,summary)
}

\keyword{iteration}