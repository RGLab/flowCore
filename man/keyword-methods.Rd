\name{keyword-methods}
\docType{methods}
\alias{keyword-methods}
\alias{keyword}
\alias{keyword,flowFrame,missing-method}

\title{Methods to retrieve keywords of a flowFrame}
\description{
  Accessor methods for items in the description slot (usually the items
  of a FCS file header). It lists the \code{keyword} and its  values for a flowFrame.

}
\usage{
keyword(object, keyword)
}
\arguments{
  \item{object}{Object of class \code{\link[flowCore]{flowFrame}}.}
  \item{keyword}{character vector or list of potential keywords. If missing the methods return all the keywords}
}
\author{N. LeMeur}
\examples{
 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 keyword(samp)

 keyword(samp,"FCSversion")

}
\keyword{methods}

