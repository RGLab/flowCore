% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/identifier-methods.R
\docType{methods}
\name{identifier-methods}
\alias{identifier-methods}
\alias{identifier}
\alias{identifier,filter-method}
\alias{identifier,transform-method}
\alias{identifier,filterReference-method}
\alias{identifier,NULL-method}
\alias{identifier,filterResult-method}
\alias{identifier,flowFrame-method}
\alias{identifier<-}
\alias{identifier<-,flowFrame-method}
\alias{identifier<-,flowFrame,ANY-method}
\title{Retrieve the GUID of flowCore objects}
\usage{
identifier(object)
}
\arguments{
\item{object}{Object of class \code{\linkS4class{flowFrame}},
\code{\linkS4class{filter}} or \code{\linkS4class{filterResult}}.}
}
\value{
Character vector representing the GUID or the name of the file.
}
\description{
Retrieve the GUID (globally unique identifier) of a \code{\link{flowFrame}}
that was generated by the cytometer or the identifier of a
\code{\linkS4class{filter}} or \code{\linkS4class{filterResult}} given by
the analyst.
}
\details{
GUID or Globally Unique Identifier is a pseudo-random number used in
software applications. While each generated GUID is not guaranteed to be
unique, the total number of unique keys (2\^128) is so large that the
probability of the same number being generated twice is very small.

Note that if no GUID has been recorded along with the FCS file, the name of
the file is returned.
}
\section{Methods}{


\describe{

\item{identifier(object = "filter")}{Return identifier of a \code{\linkS4class{filter}}
object.}

\item{identifier(object = "filterReference")}{Return identifier of a
\code{\linkS4class{filterReference}} object.}

\item{identifier(object = "filterResult")}{Return identifier of a
\code{\linkS4class{filterResult}} object.}

\item{identifier(object = "transform")}{Return identifier of a
\code{\linkS4class{transform}} object.}

\item{identifier(object = "flowFrame")}{Return GUID from the \code{description} slot of
a \code{\linkS4class{flowFrame}} object or, alternatively, the name of the
input \acronym{FCS} file in case none can be found. For
\code{\linkS4class{flowFrame}} objects there also exists a replacement
method.}

}
}

\examples{

 samp <- read.FCS(system.file("extdata","0877408774.B08", package="flowCore"))
 identifier(samp)


}
\author{
N. LeMeur
}
\keyword{methods}
