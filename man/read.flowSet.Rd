\name{read.flowSet}
\alias{read.flowSet}
\title{Read a set of FCS files}
\description{Read one or several FCS files: Data File Standard for Flow
  Cytometry}
\usage{
   read.flowSet(files=NULL,path=".",pattern=NULL,phenoData,descriptions,sep="\t",\dots)
}

\arguments{
  \item{files}{Optional character vector with filenames.}
  \item{path}{Directory where to look for the files.}
  \item{pattern}{This argument is passed on to \code{\link[base]{dir}} (see
    details). }
  \item{phenoData}{Either an object of class \code{phenoData} or 
    \code{character}.}
   \item{descriptions}{Character vector to annotate the object of class \code{\link{flowSet-class}{flowSet}}.}
  \item{sep}{Separator character that gets passed on to
    \code{\link[Biobase::read.phenoData]{read.phenoData}}.}
  \item{\dots}{Further arguments that get passed on to
    \code{\link[Biobase::read.phenoData]{read.phenoData}}, see details.}
}

\details{There are three different ways to specify the file names:
  
  First, if the argument \code{phenoData} is present and is of class
  \code{\link[Biobase:phenoData]{phenoData}}, then it is
  obtained from its column \code{name}. The column is mandatory, and an
  error will be generated if it is not there.
  Alternatively, the argument \code{phenoData} can be of class
  \code{character}, in which case this function tries to read a
  \code{phenoData} object from the file with that name by 
  calling
  \code{\link[Biobase::read.phenoData]{read.phenoData(file.path(path,
      phenoData),\dots)}}.

  Second, if the argument \code{phenoData} is not present and
  the argument \code{files} is not \code{NULL}, then \code{files}
  is expected to be a character vector with the file names.
  
  Third, if neither the argument \code{phenoData} is present nor
  \code{files} is not \code{NULL}, then the file names are obtained by
  calling   \code{dir(path, pattern)}.
}

\value{An object of class \code{\link{flowSet-class}{flowSet}}.}

\author{F. Hahne, N.Le Meur}

\examples{

fcs.loc <- system.file("extdata",package="flowCore")
file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")

samp <- read.flowSet(file.location[1:3])

}

\keyword{IO}

