\name{spillover_ng}
\alias{spillover_ng}
\alias{spillover_ng,flowSet-method}

\title{Compute a spillover matrix from a flowSet, simplified API}
\description{
	Spillover information for a particular experiment is often 
obtained by running several tubes of beads or cells stained with
a single color that can then be used to determine a spillover
matrix for use with \code{\link{compensate}}.

Matching stain channels to compensation controls is done via a matching csv file with columns 'filename'
and 'channel'. The 'channel' entries should exactly match the channel names in the FCS file. The 'filename' should be the FCS file name of each compensation control as well as one unstained control with the channel 'unstained'. 

Pregating is always done on the channels using this API, and the mode of the channel is used to compute
the spillover matrix. FSC and SSC channels can be provided to allow a pregating on (approximately) a population in the FSC and SSC dimensions.
}
\usage{
\S4method{spillover_ng}{flowSet}(x, fsc = "FSC-A", 
        ssc = "SSC-A",
       plot = FALSE, matchfile = NULL, \dots)
}

\arguments{
  
  \item{x}{A flowSet of compensation beads or cells }
  
  \item{fsc}{The name or index of the forward scatter parameter}
  
  \item{ssc}{The name or index of the side scatter parameter}
  
  \item{plot}{logical. Plots the kernel density for each channel when
  pregating. Displays the gate used. If \code{pregate} is set to \code{FALSE},
  this argument is ignored.}

  \item{matchfile}{Name of the csv file holding the compensation control file to channel matching infomration.}
    
    \item{\dots}{Additional arguments passed to \code{\link[flowStats]{rangeGate}}.}

}

\details{
  
  The algorithm used is fairly simple. First, using the scatter
  parameters, we restrict ourselves to the most closely clustered
  population to reduce the amount of debris. The selected statistic is
  then calculated on all appropriate parameters and the unstained values
  swept out of the matrix.  Every sample is then normalized to [0,1]
  with respect to the maximum value of the sample, giving the spillover
  in terms of a proportion of the primary channel intensity.

}

\value{A matrix for each of the parameters}
\references{C. B. Bagwell \& E. G. Adams (1993). Fluorescence spectral overlap compensation 
for any number of flow cytometry parameters. in: Annals of the New York Academy of Sciences, 
677:167-184.}
\author{B. Ellis}
\seealso{\code{\link{compensate}}, \code{\link{spillover}}}


\keyword{methods}  
