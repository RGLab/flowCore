\name{getIndexSort}

\docType{methods}
\alias{getIndexSort}
\alias{getIndexSort,flowFrame-method}

\title{Extract Index Sorted Data from an FCS File}


\description{

  Retrieve a data frame of index sorted data and sort indices from an FCS file.

}

\section{Methods}{

  \describe{

    \item{x = "flowFrame"}{Return a matrix of
    fluorescence intensities and indices into the sorting plate for each cell.}

  }
}


\details{
  The input FCS file should already be compensated.
  Index sorting permits association of cell-level fluorescence intensities with 
  downstream data collection on the sorted cells. Cells are sorted into a plate with \code{X,Y} 
  coordinates, and those coordinates are stored in the FCS file.
  
  This function will extract the data frame of flow data and the \code{X,Y} coordinates for the
  cell-level data, which can be written to a text file, or concatenated with sample-level 
  information and analyzed in R. The coordinates are names 'XLoc','YLoc', and a 'name' column is
  also prepended with the FCS file name.

}

\value{

  Matrix of fluorescence intensities and sort indices for plate location.
  When no index sorting data is available, invisibly returns 0. Test for 0 to check success.

}


\author{G. Finak}

\examples{

 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 # This will return a message that no index sorting data is available
 getIndexSort(samp)

}
\keyword{methods}  


