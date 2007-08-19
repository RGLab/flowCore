\name{read.FCS}
\alias{read.FCS}
\title{Read an FCS file}
\description{Read Data File Standard for Flow Cytometry}
\usage{
   read.FCS(filename, transformation = "linearize", which.lines=NULL,
                      debug = FALSE, alter.names = FALSE, column.pattern = NULL)
}

\arguments{
  \item{filename}{Character of length 1: filename}
  \item{transformation}{An character string that defines the
    type of transformation. Valid values are \code{linearize} (default) or
    \code{scale}.The  \code{linearize} tramsformation applies the appropriate power
    transform to the data while the code{scale} transformation scales all columns to [0,1]
    as in the FCS4 specification. A logical can also be used: \code{TRUE} is
    equal to \code{linearize} and \code{FALSE} correponds to no
    transformation.}
  \item{which.lines}{Numeric vector to specify the indices of the
    lines to be read. If NULL all the records are read.}
  \item{debug}{boolean indicating whether or not to print the
    debugging statements, default is TRUE}
  \item{alter.names}{boolean indicating whether or not we should rename
    the columns to valid R names using \code{\link{make.names}}. The default is FALSE.}
  \item{column.pattern}{An optional regular expression defining parameters we should keep when loading the file. The default is NULL. }
}

\details{
The function \code{read.FCS} works with the output of the FACS machine
software from a number of vendors (FCS 2.0, FCS 3.0 and List Mode Data
LMD). However, the FCS 3.0 standard includes some options that
are not yet implemented in this function. If you need extensions,
please let me know. The output of the function is an object of class
\code{flowFrame}.

For specifications of FCS 3.0 see
\url{http://www.isac-net.org} and the file
\url{../doc/fcs3.html} in the \code{doc} directory of the package.

The \code{nlines} and \code{sampling} arguments allow you to read a
  subset of the record as you might not want to read the thousands of
  events recorded in the FCS file.

The \code{which.lines} argument allows you to read a specific number of
records. 
  
}

\value{An object of class \code{\link[flowCore]{flowFrame}} that
  contains the data in the \code{exprs} slot, the parameters monitored
  in the \code{parameters} slot and the keywords and value saved in the
  header of the FCS file.}

\author{F. Hahne, N.Le Meur}

\seealso{\code{link[flowCore]{read.flowSet}}}

\examples{
samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"), transformation="linearize")
exprs(samp[1:3,])
description(samp)[3:6]
class(samp)

subset <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"), which.lines=c(2:5), transformation="linearize")
subset@exprs

theLines <- sample(1:800, 10)
subset2 <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"), which.lines=theLines, transformation="linearize")
subset2@exprs

}

\keyword{IO}

