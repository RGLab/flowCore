\name{FCSTransTransform}
\alias{FCSTransTransform}
\title{Computes a transform using the 'iplogicle' function}
\description{
	Transforms FCS data using the iplogicle function from FCSTrans by Quian et al.
  The core functionality of FCSTrans has been imported to produce transformed
  FCS data rescaled and truncated as produced by FCSTrans. The \code{w}
  parameter is estimated by \code{iplogicle} automatically, then makes a call to
  \code{iplogicore} which in turn uses the logicle transform code of Wayne
  Moore.

  For the details of the FCSTrans transformation, we recommend the excellent
  Supplementary File that accompanies Quian et al. (2012):
  \url{http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22037/suppinfo}
}
\usage{
FCSTransTransform(transformationId = "defaultFCSTransTransform", 
			channelrange = 2^18, channeldecade = 4.5, 
			range = 4096, cutoff = -111, w = NULL, rescale = TRUE)
}

\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used
  by the transform/filter  routines. } 
  \item{channelrange}{ is the range of the data. By default, 2^18 = 262144.}
  \item{channeldecade}{ is the number of logarithmic decades. By default,
    it is set to 4.5.}
  \item{range}{ the target resolution. The default value is 2^12 = 4096.}
  \item{cutoff}{ a threshold below which the logicle transformation maps values
  to 0.}
  \item{w}{ the logicle width. This is estimated by \code{iplogicle} by
  default. Details can be found in the Supplementary File from Quian et
  al.}
  \item{rescale}{ logical parameter whether or not the data should be
    rescaled to the number of channels specified in \code{range}. By
    default, the value is \code{TRUE} but can be set to FALSE
    if you want to work on the transformed scale.}
}

\references{Y Quian, Y Liu, J Campbell, E Thompson, YM Kong, RH Scheuermann; FCSTrans: An open source software system for FCS file conversion and data transformation. Cytometry A, 2012 }
\author{ Wayne Moore, N Gopalakrishnan}

\seealso{\code{\link[flowCore]{inverseLogicleTransform}}, 
    \code{\link[flowCore]{estimateLogicle} },
     \code{\link[flowCore]{logicleTransform}}}
\examples{
data(GvHD)
samp <- GvHD[[1]] 
## User defined logicle function
lgcl <- transformList(c('FL1-H', 'FL2-H'), FCSTransTransform())
after <- transform(samp, lgcl)
}

\keyword{methods}
