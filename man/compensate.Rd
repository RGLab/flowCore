\name{compensate}
\alias{compensate}
\alias{compensate,flowFrame-method}
\alias{compensate,flowSet-method}
\alias{compensate,matrix-method}


\title{Apply a compensation matrix  on a flowSet or a flowFrame object}
\description{Apply a compensation matrix  on a flowSet or a flowFrame
  object to take into account visualization. 
}

\usage{
compensate(x,spillover)
}
\arguments{
  \item{x}{Object of class \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.}
   \item{spillover}{matrix used to store spillover between
     channels. If this matrix has been store in a
    }
   }

\details{
  Note that for a  \code{\link[flowCore]{flowFrame}} object, the
  spillover matrix might have been stored under the  SPILL keyword
  and can be retrieved by using the \code{\link[flowCore]{spillover}} method.
  This matrix can also be calculated from a \code{\link[flowCore]{flowSet}} using an equivalen \code{\link[flowCore]{spillover} method.
    Finally, It can also come from a third party software.}
}


\value{
  Object of class \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}.
}
  

\seealso{
  \code{\link[flowCore]{spillover}}
}

\author{B. Ellis, N. Le Meur}
\examples{
\dontrun{
  samp <- read.FCS(system.file("extdata",
   "0877408774.B08", package="flowCore"))

 fcs.loc <- system.file("extdata/compdata", package="flowCore")
 file.location <- paste(fcs.loc, dir(fcs.loc), sep="/")
 comp.mat <-read.table(file.location[1],header=TRUE,skip=2)

 comp.samp <- compensate(samp,comp.mat)

 }
}
\keyword{methods}  

