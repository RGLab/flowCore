\name{read.FCSheader}
\alias{read.FCSheader}
\title{Read the TEXT section of a FCS file}
\description{Read (part of) the TEXT section of a Data File Standard for Flow
  Cytometry that contains FACS keywords.}
\usage{
   read.FCSheader(files, path=".", keyword=NULL, emptyValue = TRUE)
}

\arguments{
  \item{files}{Character vector of filenames.}
  \item{path}{Directory where to look for the files.}
  \item{keyword}{An optional character vector that specifies the FCS
    keyword to read.}
  \item{emptyValue}{see \code{link[flowCore]{read.FCS}}}
}

\details{
The function \code{read.FCSheader} works with the output of the FACS machine
software from a number of vendors (FCS 2.0, FCS 3.0 and List Mode Data
LMD). The output of the function is the TEXT section of the
FCS files. The user can specify some keywords to limit the output to the 
information of interest.
}

\value{A list of character vector. Each element of the list correspond
  to one FCS file.}

\author{N.Le Meur}

\seealso{\code{link[flowCore]{read.flowSet}}, \code{link[flowCore]{read.FCS}}}

\examples{
samp <- read.FCSheader(system.file("extdata",
   "0877408774.B08", package="flowCore"))
samp

samp <- read.FCSheader(system.file("extdata",
   "0877408774.B08", package="flowCore"), keyword=c("$DATE", "$FIL"))
samp
}

\keyword{IO}

