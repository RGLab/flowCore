\name{read.FCS}
\alias{read.FCS}
\title{Read an FCS file}
\description{Read Data File Standard for Flow Cytometry}
\usage{
   read.FCS(filename, transformation = "linearize", debug = FALSE, alter.names=FALSE, column.pattern=NULL)
}

\arguments{
  \item{filename}{Character of length 1: filename}
  \item{transformation}{An character string that defines the
    type of transformation. Valid values are \code{linearize} (default) or
    \code{standardize}.The former applies the appropriate power
    transform to the data while the latter scales all columns to [0,1]
    as in the FCS4 specification. A logical can also be used: \code{TRUE} is
    equal to \code{linearize} and \code{FALSE} correponds to no transformation.}
  \item{debug}{boolean indicating whether or not to print the
    debugging statements, default is TRUE}
  \item{alter.names}{boolean indicating whether or not we should rename the columns to valid R names using \code{\link{make.names}}. The default is FALSE.}
  \item{column.pattern}{An optional regular expression defining parameters we should keep when loading the file. The default is NULL. }
}

\details{
The function \code{read.FCS} works with the output of the FACS machine
software from a number of vendors. However, the FCS 3.0 standard
includes some options that are not yet implemented in this function. If
you need extensions, please let me know. The output of the function is
an object of class \code{flowFrame}.


For specifications of FCS 3.0 see
\url{http://www.isac-net.org} and the file
\url{../doc/fcs3.html} in the \code{doc} directory of the package.
}

\value{An object of class \code{\link[flowCore]{flowFrame}}.}

\author{F. Hahne, N.Le Meur}


\examples{
samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"), transformation="linearize")
exprs(samp[1:3,])
description(samp)[3:6]
class(samp)
}

\keyword{IO}

