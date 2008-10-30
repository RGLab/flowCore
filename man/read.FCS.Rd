\name{read.FCS}
\alias{read.FCS}
\alias{cleanup}
\title{Read an FCS file}

\description{Read Data File Standard for Flow Cytometry}

\usage{

   read.FCS(filename, transformation="linearize", which.lines=NULL,
            debug=FALSE, alter.names=FALSE, column.pattern=NULL,
            decades=0, ncdf = FALSE)

   cleanup()

}

\arguments{
  
  \item{filename}{Character of length 1: filename}
  
  \item{transformation}{An character string that defines the type of
    transformation. Valid values are \code{linearize} (default) or
    \code{scale}.The \code{linearize} tramsformation applies the
    appropriate power transform to the data while the code{scale}
    transformation scales all columns to $[0,10^decades]$.  defaulting
    to decades=0 as in the FCS4 specification. A logical can also be
    used: \code{TRUE} is equal to \code{linearize} and \code{FALSE}
    correponds to no transformation.}
  
  \item{which.lines}{Numeric vector to specify the indices of the lines
    to be read. If NULL all the records are read, if of length 1, a
    random sample of the size indicated by \code{which.lines} is read
    in.}

  \item{debug}{boolean indicating whether or not to print the
    debugging statements, default is TRUE}
  
  \item{alter.names}{boolean indicating whether or not we should rename
    the columns to valid R names using \code{\link{make.names}}. The
    default is FALSE.}
  
  \item{column.pattern}{An optional regular expression defining
  parameters we should keep when loading the file. The default is
  NULL. }

  \item{decades}{When scaling is activated, the number of decades to use
    for the output.}

  \item{ncdf}{Instead of reading all data into memory, this switches to
    file-based data storage. A \code{netCDF} file is creates for each
    \code{\link{flowFrame}} in the \code{.flowCoreNcdf}
    subdirectory. For large data sets this significantly reduces the
    memory profile of the R session, to the cost of speed and disk
    space. The \code{\link{exprs}} and \code{\link{exprs<-}} methods
    make sure that the user always gets a matrix of data values.  Please
    note that currently all operations that call \code{\link{exprs<-}},
    either explicitely or implicitely, will result in the creation of a
    new \code{netCDF} file. This behaviour may change in the
    future. Currently the software does not remove any of the
    \code{netCDF} files and it is up to the user to do clean up. The
    easiest way to do that is to delete the whole \code{netCDF}
    directory. To this end, one can envoke the \code{cleanup}
    function. }
  
}

\details{
  
  The function \code{read.FCS} works with the output of the FACS machine
  software from a number of vendors (FCS 2.0, FCS 3.0 and List Mode Data
  LMD). However, the FCS 3.0 standard includes some options that are not
  yet implemented in this function. If you need extensions, please let
  me know. The output of the function is an object of class
  \code{flowFrame}.

  For specifications of FCS 3.0 see \url{http://www.isac-net.org} and
  the file \url{../doc/fcs3.html} in the \code{doc} directory of the
  package.

  The \code{nlines} and \code{sampling} arguments allow you to read a
  subset of the record as you might not want to read the thousands of
  events recorded in the FCS file.

  The \code{which.lines} argument allows you to read a specific number
  of records.
  
}

\value{

  An object of class \code{\link[flowCore]{flowFrame}} that contains the
  data in the \code{exprs} slot, the parameters monitored in the
  \code{parameters} slot and the keywords and value saved in the header
  of the FCS file.

}

\author{F. Hahne, N.Le Meur}

\seealso{\code{link[flowCore]{read.flowSet}}}

\examples{
## a sample file
fcsFile <- system.file("extdata", "0877408774.B08", package="flowCore")

## read file and linearize values
samp <-  read.FCS(fcsFile, transformation="linearize")
exprs(samp[1:3,])
description(samp)[3:6]
class(samp)

## Only read in lines 2 to 5
subset <- read.FCS(fcsFile, which.lines=2:5, transformation="linearize")
exprs(subset)

## Read in a random sample of 100 lines
subset <- read.FCS(fcsFile, which.lines=100, transformation="linearize")
nrow(subset)

}

\keyword{IO}

