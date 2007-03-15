\name{spillover}
\alias{spillover}{flowSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a spillover matrix from a flowSet}
\description{
	Spillover information for a particular experiment is often 
obtained by running several tubes of beads or cells stained with
a single color that can then be used to determine a spillover
matrix for use with \code{\link{compensate}}. 
}
\usage{
spillover(x, unstained = NULL, patt = NULL, fsc = "FSC-A", 
        ssc = "SSC-A", method = "median")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A flowSet of compensation beads or cells }
  \item{unstained}{The name of index of the unstained negative 
control}
  \item{patt}{An optional regular expression defining which parameters 
should be considered.}
  \item{fsc}{The name or index of the forward scatter parameter}
  \item{ssc}{The name or index of the side scatter parameter}
  \item{method}{The statistic to use for calculation. Traditionally,
this has been the median so it is the default. The mean is sometimes
more stable.}
}
\details{
The algorithm used is fairly simple. First, using the scatter parameters,
we restrict ourselves to the most closely clustered population to reduce
the amount of debris. The selected statistic is then calculated on all 
appropriate parameters and the unstained values swept out of the matrix.
Every sample is then normalized to [0,1] with respect to the maximum value
of the sample, giving the spillover in terms of a proportion of the primary
channel intensity.
}
\value{A matrix for each of the parameters}
\references{ ~put references to the literature/web site here ~ }
\author{B. Ellis}
\seealso{\code{\link{compensate}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x,...) standardGeneric("spillover")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}  
