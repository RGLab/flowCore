\name{FCSTransTransform}
\alias{FCSTransTransform}
\title{Computes a transform using the 'iplogicle' function}
\description{
	Transforms FCS data using the iplogicle function from FCSTrans by Quian et. al. The core functionality of FCSTrans has been imported to produce transformed FCS data rescaled and truncated as produced by FCSTrans. The w parameter is estimated by 'iplogicle' automatically, then makes a call to 'iplogicore' which in turn uses the logicle transform code of Wayne Moore.
}
\usage{
FCSTransTransform(transformationId="defaultFCSTransTransform",channelrange=16777215,channeldecade=7.224719870049579 ,rescale=TRUE)
}

\arguments{
  \item{transformationId}{ A name to assign to the transformation. Used
  by the transform/filter  routines. } 
  \item{channelrange}{ is the range of the data. By default 16777215, as in FCSTrans.}
  \item{channeldecade}{ is the number of logarithmic decades. By default it is set to the value in FCSTrans 7.224719870049579.}
  \item{rescale}{ logical parameter whether or not the data should be rescaled to 4096 channels. TRUE by default, but can be set to FALSE if you want to work on the transformed scale. 
}

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
lgcl <- FCSTransTransform()
after <- transform(samp, `FL1-H` = lgcl(`FL1-H`), `FL2-H` = lgcl(`FL2-H`))
}

\keyword{methods}
