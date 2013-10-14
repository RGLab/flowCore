\name{write.flowSet}
\alias{write.flowSet}
\title{Write an FCS file}

\description{Write FCS file for each flowFrame in a flowSet}

\usage{
   write.flowSet(x, outdir=identifier(x), filename, \dots)
}

\arguments{
  
  \item{x}{ A \code{\link[flowCore:flowSet-class]{flowSet}}.}

  \item{filename}{A character scalar or vector giving the output file names. By
    default, the function will use the identifiers of the individual \code{flowFrames} as
    the file name, potentially adding the \code{.fcs} suffix unless a file
    extension is already present. Alternatively, one can supply either a
    character scalar, in which case the prefix \code{i_} is appended
    (\code{i} being an integer in \code{seq_len(length(x))}), or a
    character vector of the same length as the \code{flowSet x}.}
  
  \item{outdir}{A character scalar giving the output directory. As the
    default, the output of \code{identifier(x)} is used.  }

  \item{\dots}{Further arguments that are passed on to
    \code{\link[flowCore]{write.FCS}}.}
  
}

\details{
  
  The function \code{write.flowSet} creates FCS 3.0 standard file for
  all \code{flowFrames} in an object of class \code{flowSet}. In
  addition, it will write the content of the \code{phenoData} slot in
  the ASCII file \code{"annotation.txt"}. This file can subsequently be
  used to reconstruct the whole \code{flowSet} using the
  \code{\link[flowCore]{read.flowSet}} function, e.g.:

  \code{read.flowSet(path=outdir, phenoData="annotation.txt"}

  The function uses \code{\link[flowCore]{write.FCS}} for the actual writing of
  the FCS files.

}

\value{
  
  A character of the output directory.
  
}

\author{F. Hahne}

\seealso{\code{link[flowCore]{write.FCS}}}

\examples{
## sample data
data(GvHD)
foo <- GvHD[1:5]
outDir <- file.path(tempdir(), "foo")

## now write out into  files
write.flowSet(foo, outDir)
dir(outDir)

## and read back in
bar <- read.flowSet(path=outDir, phenoData="annotation.txt")

}

\keyword{IO}

