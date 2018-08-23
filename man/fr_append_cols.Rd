\name{fr_append_cols}
\alias{fr_append_cols}
\title{Append data columns to a flowFrame}

\description{Append data columns to a flowFrame}

\usage{
   fr_append_cols(fr, cols)
}

\arguments{
  
  \item{fr}{ A \code{\link[flowCore:flowFrame-class]{flowFrame}}.}

  \item{cols}{A numeric matrix containing the new data columns to be added.
              Must has column names to be used as new channel names.}
}

\details{
  It is used to add extra data columns to the existing flowFrame. 
  It handles keywords and parameters properly to ensure the new flowFrame can be written as a valid FCS through   the function \code{write.FCS} .


  
}

\value{
  
  A flowFrame
  
}

\author{Mike Jiang}


\examples{
  data(GvHD)
  tmp <- GvHD[[1]]
  
  kf <- kmeansFilter("FSC-H"=c("Pop1","Pop2","Pop3"), filterId="myKmFilter")
  fres <- filter(tmp, kf)
  cols <- as.integer(fres@subSet)
  cols <- matrix(cols, dimnames = list(NULL, "km"))
  tmp <- fr_append_cols(tmp, cols)
  
  tmpfile <- tempfile()
  write.FCS(tmp, tmpfile) 

}

\keyword{IO}

