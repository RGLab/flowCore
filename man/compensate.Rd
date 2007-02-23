\name{compensate}
\alias{compensate}
\alias{compensate,flowFrame-method}
\alias{compensate,flowSet-method}
\alias{compensate,matrix-method}

\title{Compensate for spillover between channels in a flowSet or a flowFrame}
\description{Apply a spillover matrix to a \code{flowSet} or a \code{flowFrame} assuming a simple linear combination of values.}

\usage{
compensate(x,spillover)
}
\arguments{
  \item{x}{Object of class \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.}
   \item{spillover}{matrix used to store spillover between channels.}
}
\details{
The essential premise of compensation is that some fluorochromes may register signals in detectors that
do not correspond to their primary detector (usually a photomultiplier tube). To compensate for this fact,
some sort of standard is used to obtain the background signal (no dye) and the amount of signal on secondary
channels for each fluorochrome relative to the signal on their primary channel. 

To calculate the spillover percentage we use either the mean or the median (more often the latter) of the 
secondary signal minus the background signal for each dye to obtain \code{n} by \code{n} matrix, \code{S}, of
so-called spillover values, expressed as a percentage of the primary channel. The observed values are then 
considered to be a linear combination of the true fluorescence and the spillover from each other channel so we
can obtain the true values by simply multiplying by the inverse of the spillover matrix.

The spillover matrix can be obtained through several means. Some flow cytometers provide a spillover matrix
calculated during acquisition, possibly by the operator, that is made available in the metadata of the flowFrame. 
While there is a theoretical standard keyword \code{\$SPILL} it can also be found in the \code{SPILLOVER} or
\code{SPILL} keyword depending on the cytometry. More commonly the spillover matrix is calculated using a series
of compensation cells or beads collected before the experiment. If you have set of FCS files with one file per 
fluorochrome as well as an unstained FCS file you can use the \code{\link[flowCore]{spillover}} method for 
\code{\link{flowSet}} to automatically calculate a spillover matrix.

}
\value{
  Object of class \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.
}

\seealso{
  \code{\link[flowCore]{spillover}}
}

\author{B. Ellis, N. Le Meur}
\examples{
 samp   <- read.flowSet(path=system.file("extdata","compdata","data",package="flowCore"))
 cfile <- system.file("extdata","compdata","compmatrix", package="flowCore")
 comp.mat <- read.table(cfile, header=TRUE, skip=2, check.names = FALSE)
 comp.mat
 summary(samp)
 samp <- compensate(samp, as.matrix(comp.mat))
 summary(samp)
}
\keyword{methods}  

