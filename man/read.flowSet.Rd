\name{read.flowSet}
\alias{read.flowSet}
\title{Read a set of FCS files}
\description{Read one or several FCS files: Data File Standard for Flow
  Cytometry}
\usage{
   read.flowSet(files=NULL, path=".",pattern=NULL,phenoData,descriptions,name.keyword,sep="\t",\dots)
}

\arguments{
  \item{files}{Optional character vector with filenames.}
  \item{path}{Directory where to look for the files.}
  \item{pattern}{This argument is passed on to \code{\link[base]{dir}},
    see details. }
  \item{phenoData}{An object of class \code{AnnotatedDataFrame},
    \code{character} or a list of values to be extracted from the
    \code{\link[flowCore]{flowFrame}} object, see
    details.}
  \item{descriptions}{Character vector to annotate the object of class
    \code{\link{flowSet-class}{flowSet}}.}
  \item{name.keyword}{An optional character vector that specifies which
    FCS keyword to use as the sample names.}
  \item{sep}{Separator character that gets passed on to
    \code{\link[Biobase::read.AnnotatedDataFrame]{read.AnnotatedDataFrame}}.}
  \item{\dots}{Further arguments that get passed on to
    \code{\link[flowCore]{read.FCS}}, see details.}
}

\details{There are four different ways to specify the file from which
  data is to be imported:
  
  First, if the argument \code{phenoData} is present and is of class
  \code{\link[Biobase:AnnotatedDataFrame]{AnnotatedDataFrame}}, then the
  file names are obtained from its sample names (i.e. column \code{name}. The column is mandatory, and an
  error will be generated if it is not there.
  Alternatively, the argument \code{phenoData} can be of class
  \code{character}, in which case this function tries to read a
  \code{AnnotatedDataFrame} object from the file with that name by 
  calling
  \code{\link[Biobase::read.AnnotatedDataFrame]{read.AnnotatedDataFrame(file.path(path,
      phenoData),\dots)}}.
  
  In some cases the file names are not a reasonable selection criterion
  and the user might want to import files based on some keywords within
  the file. One or several keyword value pairs can be given as the phenoData
  argument in form of a named list. 

  Third, if the argument \code{phenoData} is not present and
  the argument \code{files} is not \code{NULL}, then \code{files}
  is expected to be a character vector with the file names.
  
  Third, if neither the argument \code{phenoData} is present nor
  \code{files} is not \code{NULL}, then the file names are obtained by
  calling   \code{dir(path, pattern)}.
}

\value{An object of class \code{\link{flowSet-class}{flowSet}}.}

\author{F. Hahne, N.Le Meur, B. Ellis}

\examples{

fcs.loc <- system.file("extdata",package="flowCore")
file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")

samp <- read.flowSet(file.location[1:3])

}

\keyword{IO}

