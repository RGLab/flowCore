\name{read.flowSet}
\alias{read.flowSet}
\title{Read a set of FCS files}
\description{Read one or several FCS files: Data File Standard for Flow
  Cytometry}
\usage{

read.flowSet(files=NULL, path=".", pattern=NULL, phenoData,
             descriptions,name.keyword, alter.names=FALSE,
             transformation = "linearize", which.lines=NULL,
             column.pattern = NULL, invert.pattern = FALSE, decades=0, sep="\t",
             as.is=TRUE, name, ncdf=FALSE, dataset=NULL, min.limit=NULL,
             truncate_max_range = TRUE, emptyValue=TRUE, 
             ignore.text.offset = FALSE, \dots)
             emptyValue=TRUE, ignore.text.offset = FALSE, \dots)

}

\arguments{
  
  \item{files}{Optional character vector with filenames.}
  
  \item{path}{Directory where to look for the files.}
  
  \item{pattern}{This argument is passed on to
    \code{\link[base:list.files]{dir}}, see details. }
  
  \item{phenoData}{An object of class \code{AnnotatedDataFrame},
    \code{character} or a list of values to be extracted from the
    \code{\link[flowCore:flowFrame-class]{flowFrame}} object, see
    details.}
  
  \item{descriptions}{Character vector to annotate the object of class
    \code{\link[flowCore:flowSet-class]{flowSet}}.}
  
  \item{name.keyword}{An optional character vector that specifies which
    FCS keyword to use as the sample names. If this is not set, the GUID
    of the FCS file is used for sampleNames, and if that is not present
    (or not unique), then the file names are used. }
  
  \item{alter.names}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{transformation}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{which.lines}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{column.pattern}{see \code{\link[flowCore]{read.FCS}} for
    details.}

  \item{invert.pattern}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{decades}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{sep}{Separator character that gets passed on to
    \code{\link[Biobase]{read.AnnotatedDataFrame}}.}

  \item{as.is}{Logical that gets passed on to
    \code{\link[Biobase]{read.AnnotatedDataFrame}}. This
    controls the automatic coercion of characters to factors in the
    \code{phenoData}slot. }

  \item{name}{An optional character scalar used as name of the object.}

  \item{ncdf}{Deprecated. Please refer to 'ncdfFlow' package for cdf based storage.}

  \item{dataset}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{min.limit}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{truncate.max.range}{see \code{link[flowCore]{read.FCS}} for details.}
  
  \item{emptyValue}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{ignore.text.offset}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{truncate_max_range}{see \code{\link[flowCore]{read.FCS}} for details.}
  
  \item{\dots}{Further arguments that get passed on to
    \code{\link[Biobase]{read.AnnotatedDataFrame}}, see details.}
  
}

\details{There are four different ways to specify the file from which
  data is to be imported:
  
  First, if the argument \code{phenoData} is present and is of class
  \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}, then the
  file names are obtained from its sample names (i.e. row names of the underlying data.frame).
  Also column \code{name} will be generated based on sample names if it is not there. This column is mainly used by visualization methods in flowViz.   
  Alternatively, the argument \code{phenoData} can
  be of class \code{character}, in which case this function tries to
  read a \code{AnnotatedDataFrame} object from the file with that name
  by calling
  \code{\link[Biobase]{read.AnnotatedDataFrame}(file.path(path,phenoData),\dots)}.  
  
  In some cases the file names are not a reasonable selection criterion
  and the user might want to import files based on some keywords within
  the file. One or several keyword value pairs can be given as the phenoData
  argument in form of a named list. 

  Third, if the argument \code{phenoData} is not present and
  the argument \code{files} is not \code{NULL}, then \code{files}
  is expected to be a character vector with the file names.
  
  Fourth, if neither the argument \code{phenoData} is present nor
  \code{files} is not \code{NULL}, then the file names are obtained by
  calling   \code{dir(path, pattern)}.
}

\value{An object of class \code{\link[flowCore:flowSet-class]{flowSet}}.}

\author{F. Hahne, N.Le Meur, B. Ellis}

\examples{

fcs.loc <- system.file("extdata",package="flowCore")
file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")

samp <- read.flowSet(file.location[1:3])

}

\keyword{IO}

