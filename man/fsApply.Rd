\name{fsApply}
\alias{fsApply}
\alias{fsApply,flowSet,ANY}


\title{Apply a Function over values in a flowSet}

\description{
\code{fsApply} like many of the \code{apply}-style functions in R acts as an iterator for
\code{flowSet} objects, allowing the application of a function to either the \code{flowFrame} 
or the data matrix itself. The output can the be reconstructed as either a \code{flowSet}, a list
or a matrix depending on options and the type of objects returned.
}

\usage{
    fsApply(x, FUN,\dots, simplify=TRUE, use.exprs=FALSE )
}

\arguments{
  \item{x}{\code{\link[flowCore]{flowSet}} to be used}
  \item{FUN}{the function to be applied to each element of \code{x}}
  \item{simplify}{logical (default: TRUE);  if all true and all objects are \code{flowFrame} objects, 
	a \code{flowSet} object will be constructed. If all of the values are of the same type there will be
	an attempt to construct a vector or matrix of the appropriate type (e.g. all numeric results will return a matrix).}
  \item{use.exprs}{logical (default: FALSE); should the \code{FUN} be applied on the
    \code{\link[flowCore]{flowFrame}} object or the expression values}
  \item{\dots}{optional arguments to \code{FUN}.}
  
}
\seealso{
  \code{\link{apply}}, \code{\link{sapply}}
}
\examples{
fcs.loc <- system.file("extdata",package="flowCore")
file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")
samp <- read.flowSet(file.location[1:3]) 

#Get summary information about each sample.
fsApply(samp,summary)

#Obtain the median of each parameter in each frame.
fsApply(samp,each_col,median)

}

\keyword{iteration}